import netCDF4 as nc
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import json
import logging
from datetime import datetime, timezone
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Any
import traceback
import gc as _gc


# ══════════════════════════════════════════════════════════════════════════════
# ⚠️  HIERARQUIA DE EXCEÇÕES — erros tipados por estágio do pipeline
# ══════════════════════════════════════════════════════════════════════════════


class GCAnalyzerError(Exception):
    def __init__(self, message: str, context: dict | None = None):
        super().__init__(message)
        self.context = context or {}


class CDFReadError(GCAnalyzerError):
    """Arquivo CDF ausente, corrompido ou com variáveis obrigatórias faltantes."""


class BaselineError(GCAnalyzerError):
    """Falha na subtração de baseline (Whittaker AsLS divergiu ou NaN no resultado)."""


class PeakDetectionError(GCAnalyzerError):
    """Nenhum pico sobreviveu aos filtros de SNR, largura e curvatura."""


class IntegrationError(GCAnalyzerError):
    """Falha em todas as tentativas de integração (EMG + fallback trapezoidal)."""


class AlignmentError(GCAnalyzerError):
    """Erro no alinhamento multi-corrida (IS ausente, corridas incompatíveis)."""


# ══════════════════════════════════════════════════════════════════════════════
# 📦  RunResult
# ══════════════════════════════════════════════════════════════════════════════


@dataclass
class RunResult:
    run_id: str
    status: str
    cdf_path: str
    results_df: Optional[pd.DataFrame]
    audit_events: list
    error_type: Optional[str] = None
    error_message: Optional[str] = None
    error_traceback: Optional[str] = None

    @property
    def ok(self) -> bool:
        return self.status == "OK"

    def __repr__(self) -> str:
        if self.ok:
            n = len(self.results_df) if self.results_df is not None else 0
            return f"RunResult(run_id={self.run_id!r}, status=OK, peaks={n})"
        return f"RunResult(run_id={self.run_id!r}, status=FAILED, error={self.error_type}: {self.error_message!r})"


# ══════════════════════════════════════════════════════════════════════════════
# ⚙️  PROCESSING METHOD
# ══════════════════════════════════════════════════════════════════════════════


@dataclass
class ProcessingMethod:
    name: str = "default"
    version: str = "1.0"
    description: str = ""
    created_by: str = ""
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds"))
    source_file: str = ""

    baseline_lam: float = 1e8
    baseline_p: float = 0.0001
    noise_percentile: int = 20
    snr_threshold: float = 3.0
    min_width_seconds: float = 1.0
    min_distance_seconds: float = 2.0
    rs_deconv_threshold: float = 1.2
    valley_pct_independent: float = 85.0
    valley_pct_dropline: float = 50.0
    valley_pct_skim_max: float = 25.0
    height_ratio_rider: float = 0.15
    solvent_rt_cutoff_s: float = 60.0
    solvent_area_factor: float = 5.0
    is_rt_seconds: Optional[float] = None
    is_search_window_s: float = 10.0
    rrt_bin_tolerance: float = 0.02

    # ── Fator de capacidade / seletividade ───────────────────────────────────
    # dead_time_s : tempo morto (t₀) da coluna em segundos.
    #   Necessário para calcular k' = (tR - t₀) / t₀.
    #   0.0 → recurso desativado (k' e α ficam NaN).
    #   Pode ser medido pela injeção de metano (FID) ou determinado por
    #   geometria da coluna: t₀ ≈ L / (u × (1 + k)) onde u é fluxo linear.
    dead_time_s: float = 0.0

    # ── CQI — pesos dos sub-scores ───────────────────────────────────────────
    # O Chromatographic Quality Index combina N, Rs, TF e SNR num único valor
    # 0–1. Os pesos controlam a influência relativa de cada dimensão.
    # Por padrão, todos os pesos são iguais (média geométrica ponderada).
    cqi_weight_n: float = 1.0  # peso do score de eficiência (N_plates)
    cqi_weight_rs: float = 1.0  # peso do score de resolução (Rs_usp)
    cqi_weight_tf: float = 1.0  # peso do score de simetria (Tailing Factor)
    cqi_weight_snr: float = 1.0  # peso do score de SNR

    # ── Limites de referência para normalização do CQI ───────────────────────
    cqi_n_ref: float = 5000.0  # N de referência (pontos de eficiência "perfeita")
    cqi_rs_ref: float = 1.5  # Rs mínimo aceitável (USP: ≥ 1.5 = baseline resolved)
    cqi_snr_ref: float = 10.0  # SNR de referência (sinal considerado "robusto")

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=indent)

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "ProcessingMethod":
        path = Path(path)
        data = json.loads(path.read_text(encoding="utf-8"))
        data["source_file"] = str(path.resolve())
        return cls(**data)

    @classmethod
    def default(cls) -> "ProcessingMethod":
        return cls(name="default", description="Parâmetros padrão de fábrica.")

    def __str__(self) -> str:
        src = f" | fonte: {self.source_file}" if self.source_file else ""
        return f"ProcessingMethod(name={self.name!r}, version={self.version!r}{src})"

    def __repr__(self) -> str:
        return self.__str__()


# ══════════════════════════════════════════════════════════════════════════════
# 🗒️  AUDIT TRAIL
# ══════════════════════════════════════════════════════════════════════════════


@dataclass
class AuditEvent:
    timestamp: str
    level: str
    module: str
    message: str
    original_value: Optional[Any] = None
    new_value: Optional[Any] = None
    context: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)

    def __str__(self) -> str:
        ts = self.timestamp[11:23]
        ctx = f" | ctx={self.context}" if self.context else ""
        vals = ""
        if self.original_value is not None or self.new_value is not None:
            vals = f" | {self.original_value!r} → {self.new_value!r}"
        return f"[{ts}] {self.level:<8} [{self.module}] {self.message}{vals}{ctx}"


class AuditLogger:
    def __init__(self, run_id: Optional[str] = None, echo: bool = False):
        self.run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
        self.echo = echo
        self._events: list[AuditEvent] = []
        self._py_log = logging.getLogger(f"GCAnalyzer.{self.run_id}")

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat(timespec="milliseconds")

    def _add(self, level, module, message, original_value=None, new_value=None, **ctx):
        evt = AuditEvent(
            timestamp=self._now(),
            level=level,
            module=module,
            message=message,
            original_value=original_value,
            new_value=new_value,
            context=ctx,
        )
        self._events.append(evt)
        if self.echo:
            print(evt)
        py_level = {"INFO": logging.INFO, "WARN": logging.WARNING, "ERROR": logging.ERROR, "DECISION": logging.INFO}.get(
            level, logging.DEBUG
        )
        self._py_log.log(py_level, str(evt))

    def info(self, module, message, **ctx):
        self._add("INFO", module, message, **ctx)

    def warn(self, module, message, original_value=None, new_value=None, **ctx):
        self._add("WARN", module, message, original_value=original_value, new_value=new_value, **ctx)

    def error(self, module, message, **ctx):
        self._add("ERROR", module, message, **ctx)

    def decision(self, module, message, original_value=None, new_value=None, **ctx):
        self._add("DECISION", module, message, original_value=original_value, new_value=new_value, **ctx)

    def log_method(self, method: ProcessingMethod):
        src = method.source_file if method.source_file else "<in-memory>"
        self.info(
            "Method",
            f"ProcessingMethod carregado: nome='{method.name}', versão='{method.version}', fonte='{src}'.",
            method_name=method.name,
            method_version=method.version,
            method_source=src,
            method_description=method.description,
            method_created_by=method.created_by,
            method_created_at=method.created_at,
            params=method.to_dict(),
        )

    def log_baseline(self, lam, p, mean_reduction, pts):
        self.info(
            "Baseline",
            f"Whittaker AsLS aplicado (λ={lam:.0e}, p={p}). Redução média: {mean_reduction:.1f} u.a. ({pts} pontos).",
            lam=lam,
            p=p,
            mean_reduction=mean_reduction,
            n_points=pts,
        )

    def log_noise_global(self, sigma, fallback=False):
        tag = " [FALLBACK via diff-std]" if fallback else ""
        self.info("NoiseEstimate", f"Ruído global (MAD×1.4826){tag}: σ={sigma:.4f}", sigma=sigma)

    def log_peak_rejection(self, reason, rt_s, value=None):
        self.warn(
            "PeakDetection", f"Pico rejeitado em RT={rt_s:.2f}s — motivo: {reason}.", context_rt=rt_s, rejection_value=value
        )

    def log_peaks_found(self, n, rts):
        self.info("PeakDetection", f"{n} pico(s) detectado(s) após filtragem SNR/largura/d².", n_peaks=n, retention_times_s=rts)

    def log_snr_rejection(self, n_rejected, threshold):
        if n_rejected:
            self.warn(
                "PeakDetection",
                f"{n_rejected} pico(s) rejeitado(s) por SNR local < {threshold}.",
                snr_threshold=threshold,
                n_rejected=n_rejected,
            )

    def log_overlap_decision(self, rt1, rt2, Rs, valley_pct, height_ratio, method):
        self.decision(
            "Integration",
            f"Sobreposição RT={rt1:.2f}–{rt2:.2f}s → método '{method}' selecionado "
            f"(Rs={Rs:.2f}, %Vale={valley_pct:.1f}%, h_ratio={height_ratio:.2f}).",
            original_value="UNKNOWN",
            new_value=method,
            rt1=rt1,
            rt2=rt2,
            Rs=Rs,
            valley_pct=valley_pct,
            height_ratio=height_ratio,
        )

    def log_emg_fallback(self, rt_s, area_trap, reason):
        self.decision(
            "Integration",
            f"Ajuste EMG falhou em RT={rt_s:.2f}s ({reason}). Revertido para integração trapezoidal (área={area_trap:.0f}).",
            original_value="EMG",
            new_value="TRAPEZOID",
            rt=rt_s,
            area_trap=area_trap,
            failure_reason=reason,
        )

    def log_integration(self, method, rt, area, snr, window, extra=None):
        ctx = dict(rt=rt, area=area, snr=snr, window_pts=window)
        if extra:
            ctx.update(extra)
        self.info("Integration", f"[{method}] RT={rt:.2f}s integrado: área={area:.0f}, SNR={snr:.1f}.", **ctx)

    def log_solvent_removal(self, n_before, n_after, median_area, rt_cutoff, factor):
        n_removed = n_before - n_after
        self.decision(
            "QC",
            f"Remoção de pico de solvente: {n_removed} pico(s) eliminado(s) "
            f"(RT ≤ {rt_cutoff}s ou área > {factor}× mediana={median_area:.0f}).",
            original_value=n_before,
            new_value=n_after,
            median_area=median_area,
            rt_cutoff=rt_cutoff,
            factor=factor,
        )

    def log_deconv_audit(self, area_trap, area_emg):
        pct_err = 100 * (area_emg - area_trap) / max(area_trap, 1)
        level = "WARN" if abs(pct_err) > 10 else "INFO"
        self._add(
            level,
            "Deconvolution",
            f"Área trap total={area_trap:.0f} | EMG total={area_emg:.0f} | erro={pct_err:.1f}%.",
            area_trap=area_trap,
            area_emg=area_emg,
            pct_error=pct_err,
        )

    def to_dict_list(self):
        return [e.to_dict() for e in self._events]

    def to_dataframe(self):
        rows = []
        for e in self._events:
            d = e.to_dict()
            d["context"] = json.dumps(d["context"], ensure_ascii=False)
            rows.append(d)
        return pd.DataFrame(rows)

    def to_json(self, indent=2):
        payload = {"run_id": self.run_id, "generated": self._now(), "n_events": len(self._events), "events": self.to_dict_list()}
        return json.dumps(payload, ensure_ascii=False, indent=indent, default=str)

    def summary(self):
        counts: dict[str, int] = {}
        for e in self._events:
            counts[e.level] = counts.get(e.level, 0) + 1
        return {"run_id": self.run_id, "total": len(self._events), **counts}

    def __len__(self):
        return len(self._events)

    def __repr__(self):
        s = self.summary()
        return (
            f"AuditLogger(run_id={s['run_id']!r}, total={s['total']}, "
            f"INFO={s.get('INFO',0)}, WARN={s.get('WARN',0)}, "
            f"ERROR={s.get('ERROR',0)}, DECISION={s.get('DECISION',0)})"
        )


# ══════════════════════════════════════════════════════════════════════════════
# 🔬  GCAnalyzer
# ══════════════════════════════════════════════════════════════════════════════

from scipy.signal import find_peaks as _scipy_find_peaks, peak_widths
from scipy.optimize import curve_fit
from scipy.stats import exponnorm
from scipy.integrate import trapezoid
from pybaselines import Baseline as _PyBaseline


class GCAnalyzer:
    def __init__(self, method: Optional[ProcessingMethod] = None, run_id: Optional[str] = None, echo_audit: bool = False):
        self.method = method or ProcessingMethod.default()
        self.audit = AuditLogger(run_id=run_id, echo=echo_audit)
        self.audit.log_method(self.method)

    @property
    def _m(self) -> ProcessingMethod:
        return self.method

    # ==========================================================
    # 1️⃣  LEITURA CDF
    # ==========================================================
    def read_cdf(self, cdf_filepath):
        try:
            with nc.Dataset(cdf_filepath) as gc_data:
                required_vars = ("ordinate_values", "actual_delay_time", "actual_sampling_interval")
                missing = [v for v in required_vars if v not in gc_data.variables]
                if missing:
                    raise CDFReadError(
                        f"Variáveis obrigatórias ausentes no CDF: {missing}",
                        context={"filepath": str(cdf_filepath), "missing_vars": missing},
                    )
                intensity = np.array(gc_data.variables["ordinate_values"][:])
                delay = np.array(gc_data.variables["actual_delay_time"][:])
                interval = np.array(gc_data.variables["actual_sampling_interval"][:])
                retention_time = delay + (interval * np.arange(len(intensity)))
        except CDFReadError:
            raise
        except Exception as exc:
            raise CDFReadError(f"Falha ao ler CDF '{cdf_filepath}': {exc}", context={"filepath": str(cdf_filepath)}) from exc

        self.audit.info(
            "IO",
            f"CDF lido: {cdf_filepath} | {len(intensity)} pontos | RT=[{retention_time[0]:.1f}–{retention_time[-1]:.1f}]s.",
            filepath=cdf_filepath,
            n_points=len(intensity),
            rt_start=float(retention_time[0]),
            rt_end=float(retention_time[-1]),
        )
        return retention_time, intensity

    # ==========================================================
    # 2️⃣  BASELINE
    # ==========================================================
    def remove_baseline(self, rt, intensity):
        lam, p = self._m.baseline_lam, self._m.baseline_p
        try:
            baseline_fitter = _PyBaseline(x_data=rt)
            baseline, _ = baseline_fitter.asls(intensity, lam=lam, p=p)
        except Exception as exc:
            raise BaselineError(f"Whittaker AsLS divergiu (λ={lam:.0e}, p={p}): {exc}", context={"lam": lam, "p": p}) from exc

        if not np.all(np.isfinite(baseline)):
            raise BaselineError(
                f"Baseline contém NaN/Inf após AsLS (λ={lam:.0e}, p={p}).",
                context={"lam": lam, "p": p, "n_nan": int(np.sum(~np.isfinite(baseline)))},
            )

        corrected = intensity - baseline
        corrected[corrected < 0] = 0
        self.audit.log_baseline(lam=lam, p=p, mean_reduction=float(np.mean(baseline)), pts=len(intensity))
        return corrected, baseline

    # ==========================================================
    # 3️⃣  RUÍDO GLOBAL (MAD)
    # ==========================================================
    def estimate_noise_level(self, intensity):
        percentile = self._m.noise_percentile
        threshold = np.percentile(intensity, percentile)
        baseline_region = intensity[intensity <= threshold]
        if len(baseline_region) < 10:
            baseline_region = intensity
        median = np.median(baseline_region)
        mad = np.median(np.abs(baseline_region - median))
        sigma = mad * 1.4826
        fallback = False
        if sigma <= 0:
            sigma = np.std(np.diff(intensity)) / np.sqrt(2) * 0.1
            fallback = True
        self.audit.log_noise_global(sigma=float(sigma), fallback=fallback)
        return sigma

    # ==========================================================
    # 3️⃣b RUÍDO LOCAL (por pico)
    # ==========================================================
    def estimate_local_snr(self, rt, intensity, peak_idx, left_bound, right_bound, width_factor=3.0):
        peak_width_pts = max(right_bound - left_bound, 1)
        half_window = int(width_factor * peak_width_pts)
        l_start = max(0, left_bound - half_window)
        r_end = min(len(rt) - 1, right_bound + half_window)
        left_pts = np.arange(l_start, left_bound)
        right_pts = np.arange(right_bound, r_end)
        MIN_PTS_PER_SIDE = max(10, half_window // 2)
        if len(left_pts) < MIN_PTS_PER_SIDE and len(right_pts) >= MIN_PTS_PER_SIDE:
            extra = right_pts[: MIN_PTS_PER_SIDE - len(left_pts)]
            right_pts = np.concatenate([right_pts, extra])
        elif len(right_pts) < MIN_PTS_PER_SIDE and len(left_pts) >= MIN_PTS_PER_SIDE:
            extra = left_pts[-(MIN_PTS_PER_SIDE - len(right_pts)) :]
            left_pts = np.concatenate([left_pts, extra])
        ref_idx = np.concatenate([left_pts, right_pts]).astype(int)
        if len(ref_idx) < 6:
            global_noise = self.estimate_noise_level(intensity)
            local_bl_val = (intensity[left_bound] + intensity[right_bound]) / 2.0
            signal = max(intensity[peak_idx] - local_bl_val, 0.0)
            self.audit.warn(
                "NoiseEstimate",
                f"Poucos pontos de referência ({len(ref_idx)}) em RT={rt[peak_idx]:.2f}s. Usando ruído global.",
                n_ref_pts=len(ref_idx),
                rt=float(rt[peak_idx]),
            )
            return signal / max(global_noise, 1e-9), local_bl_val, global_noise
        x_ref = rt[ref_idx]
        y_ref = intensity[ref_idx]
        local_noise = np.std(np.diff(y_ref)) / np.sqrt(2)
        if local_noise <= 0 or not np.isfinite(local_noise):
            local_noise = self.estimate_noise_level(intensity)
            self.audit.warn(
                "NoiseEstimate", f"Ruído local inválido em RT={rt[peak_idx]:.2f}s. Usando ruído global.", rt=float(rt[peak_idx])
            )
        coeffs = np.polyfit(x_ref, y_ref, deg=1)
        local_bl_val = np.polyval(coeffs, rt[peak_idx])
        signal = max(intensity[peak_idx] - local_bl_val, 0.0)
        snr = signal / local_noise if local_noise > 0 else 0.0
        return snr, local_bl_val, local_noise

    # ==========================================================
    # 3️⃣c INTEGRAÇÃO TRAPEZOIDAL
    # ==========================================================
    def integrate_trapezoid_segment(self, rt, intensity, left, right):
        x_seg = rt[left:right]
        y_seg = intensity[left:right]
        if len(x_seg) < 2:
            self.audit.warn(
                "Integration", f"Segmento muito curto ({len(x_seg)} pts) para integração trapezoidal.", left=left, right=right
            )
            return 0.0, x_seg, np.zeros_like(x_seg), np.zeros_like(x_seg)
        y_left, y_right = float(y_seg[0]), float(y_seg[-1])
        baseline_virtual = np.linspace(y_left, y_right, len(x_seg))
        y_above = np.maximum(y_seg - baseline_virtual, 0.0)
        area = float(trapezoid(y_above, x_seg))
        self.audit.info(
            "Integration",
            f"Integração trapezoidal: janela=[{left}:{right}], área={area:.2f}.",
            window=(left, right),
            area=area,
            y_left=y_left,
            y_right=y_right,
            max_above_bl=float(np.max(y_above)),
        )
        return area, x_seg, y_above, baseline_virtual

    # ==========================================================
    # 4️⃣  DETECÇÃO DE PICOS
    # ==========================================================
    def find_peaks(self, rt, intensity):
        snr_threshold = self._m.snr_threshold
        min_width_seconds = self._m.min_width_seconds
        min_distance_seconds = self._m.min_distance_seconds
        dt = np.mean(np.diff(rt))
        min_width_pts = max(1, int(min_width_seconds / dt))
        min_distance_pts = max(1, int(min_distance_seconds / dt))
        noise_sigma = self.estimate_noise_level(intensity)
        if noise_sigma <= 0:
            noise_sigma = np.std(intensity) * 0.01
            self.audit.warn("PeakDetection", "noise_sigma ≤ 0 após MAD — usando 1% do desvio padrão.", noise_sigma=noise_sigma)
        dynamic_prominence = snr_threshold * noise_sigma
        from scipy.ndimage import gaussian_filter1d

        smoothed = gaussian_filter1d(intensity, sigma=3)
        d1 = np.gradient(smoothed, rt)
        d2 = np.gradient(d1, rt)
        self.audit.info(
            "PeakDetection",
            f"Parâmetros de detecção: dt={dt:.4f}s, min_width={min_width_pts}pts, min_dist={min_distance_pts}pts, prominence≥{dynamic_prominence:.2f}.",
            dt=dt,
            min_width_pts=min_width_pts,
            min_distance_pts=min_distance_pts,
            prominence=dynamic_prominence,
            snr_threshold=snr_threshold,
        )
        peaks, _ = _scipy_find_peaks(intensity, prominence=dynamic_prominence, distance=min_distance_pts, width=min_width_pts)
        if len(peaks) == 0:
            msg = "Nenhum pico encontrado após find_peaks — prominence muito alta ou sinal muito ruidoso."
            self.audit.warn("PeakDetection", msg)
            raise PeakDetectionError(
                msg,
                context={
                    "snr_threshold": snr_threshold,
                    "dynamic_prominence": float(dynamic_prominence),
                    "noise_sigma": float(noise_sigma),
                },
            )
        widths, _, left_ips, right_ips = peak_widths(intensity, peaks, rel_height=0.95)
        left_base = np.maximum(0, np.floor(left_ips)).astype(int)
        right_base = np.minimum(len(intensity) - 1, np.ceil(right_ips)).astype(int)
        max_reasonable_sigma = 10
        valid = []
        for i, p in enumerate(peaks):
            if d2[p] >= 0:
                self.audit.log_peak_rejection("d²<0 não satisfeito", rt_s=float(rt[p]), value=float(d2[p]))
                continue
            sigma_est = (widths[i] * dt) / 2.355
            if sigma_est >= max_reasonable_sigma:
                self.audit.log_peak_rejection(
                    f"Largura excessiva (σ_est={sigma_est:.2f}s ≥ {max_reasonable_sigma}s)", rt_s=float(rt[p]), value=sigma_est
                )
                continue
            valid.append(i)
        peaks = peaks[valid]
        left_ips = left_ips[valid]
        right_ips = right_ips[valid]
        left_base = left_base[valid]
        right_base = right_base[valid]
        snr_values = np.array(
            [self.estimate_local_snr(rt, intensity, p, left_base[i], right_base[i])[0] for i, p in enumerate(peaks)]
        )
        n_rejected = int((snr_values < snr_threshold).sum())
        self.audit.log_snr_rejection(n_rejected, snr_threshold)
        mask = snr_values >= snr_threshold
        peaks = peaks[mask]
        left_ips = left_ips[mask]
        right_ips = right_ips[mask]
        left_base = left_base[mask]
        right_base = right_base[mask]
        snr_values = snr_values[mask]
        if len(peaks) == 0:
            msg = f"Todos os picos foram rejeitados pelo filtro de SNR local < {snr_threshold}."
            self.audit.warn("PeakDetection", msg, snr_threshold=snr_threshold)
            raise PeakDetectionError(msg, context={"snr_threshold": snr_threshold})
        self.audit.log_peaks_found(n=len(peaks), rts=[round(float(rt[p]), 2) for p in peaks])
        return peaks, left_ips, right_ips, snr_values, left_base, right_base

    # ==========================================================
    # 5️⃣  MODELOS
    # ==========================================================
    @staticmethod
    def emg(x, A, mu, sigma, tau):
        if sigma <= 0 or tau <= 0:
            return np.zeros_like(x, dtype=float)
        result = A * exponnorm.pdf(x, K=tau / sigma, loc=mu, scale=sigma)
        return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)

    @staticmethod
    def multi_emg(x, *params):
        y = np.zeros_like(x, dtype=float)
        for i in range(0, len(params), 4):
            y += GCAnalyzer.emg(x, params[i], params[i + 1], params[i + 2], params[i + 3])
        return y

    @staticmethod
    def gaussian(x, A, mu, sigma):
        return A * np.exp(-((x - mu) ** 2) / (2 * sigma**2))

    @staticmethod
    def multi_gaussian(x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            y += params[i] * np.exp(-((x - params[i + 1]) ** 2) / (2 * params[i + 2] ** 2))
        return y

    # ==========================================================
    # 6️⃣  AJUSTE EMG
    # ==========================================================
    def fit_emg_peak(self, rt, intensity, peak_idx, left, right):
        x = rt[left:right]
        y = intensity[left:right]
        if len(x) < 5:
            self.audit.warn(
                "Integration",
                f"Janela muito pequena ({len(x)} pts) em RT={rt[peak_idx]:.2f}s.",
                rt=float(rt[peak_idx]),
                window_size=len(x),
            )
            return None
        area_trap, _, y_above, _ = self.integrate_trapezoid_segment(rt, intensity, left, right)
        height_above_bl = float(np.max(y_above)) if len(y_above) > 0 else 0.0
        A0 = trapezoid(y, x)
        mu0 = rt[peak_idx]
        sigma0 = max((rt[right] - rt[left]) / 6, 0.01)
        tau0 = sigma0
        emg_params = {}
        try:
            popt, _ = curve_fit(
                self.emg,
                x,
                y,
                p0=[A0, mu0, sigma0, tau0],
                bounds=([0, min(x), 0.001, 0.001], [np.inf, max(x), np.inf, np.inf]),
                maxfev=5000,
            )
            A, mu, sigma, tau = popt
            sigma, tau = abs(sigma), abs(tau)
            y_emg = self.emg(x, A, mu, sigma, tau)
            area_emg = float(trapezoid(y_emg, x))
            self.audit.log_integration(
                method="EMG",
                rt=float(mu),
                area=area_trap,
                snr=0.0,
                window=(left, right),
                extra=dict(area_emg=area_emg, sigma=sigma, tau=tau, tailing=tau / sigma if sigma > 0 else None),
            )
            emg_params = {
                "A_param": A,
                "rt": mu,
                "height": height_above_bl,
                "area": area_trap,
                "area_emg": area_emg,
                "sigma": sigma,
                "tau": tau,
                "tailing_factor": tau / sigma if sigma > 0 else np.nan,
                "tau_bounded": abs(tau - 0.001) < 1e-6,
                "marker_rt": rt[peak_idx],
                "marker_height": intensity[peak_idx],
                "peak_index_apex": int(peak_idx),
                "peak_index_start": int(left),
                "peak_index_end": int(right),
            }
        except Exception as e:
            self.audit.log_emg_fallback(rt_s=float(rt[peak_idx]), area_trap=area_trap, reason=str(e))
            emg_params = {
                "A_param": np.nan,
                "rt": float(rt[peak_idx]),
                "height": height_above_bl,
                "area": area_trap,
                "area_emg": np.nan,
                "sigma": np.nan,
                "tau": np.nan,
                "tailing_factor": np.nan,
                "tau_bounded": np.nan,
                "marker_rt": rt[peak_idx],
                "marker_height": intensity[peak_idx],
                "peak_index_apex": int(peak_idx),
                "peak_index_start": int(left),
                "peak_index_end": int(right),
            }
        try:
            popt_g, _ = curve_fit(
                self.gaussian,
                x,
                y,
                p0=[np.max(y), mu0, sigma0],
                bounds=([0, min(x), 0.001], [np.inf, max(x), np.inf]),
                maxfev=5000,
            )
            Ag, mug, sigmag = popt_g
            emg_params.update(
                {
                    "gauss_A": Ag,
                    "gauss_mu": mug,
                    "gauss_sigma": abs(sigmag),
                    "area_gauss": float(trapezoid(self.gaussian(x, Ag, mug, abs(sigmag)), x)),
                }
            )
        except Exception as e:
            self.audit.warn(
                "Integration", f"Ajuste Gaussiano falhou em RT={rt[peak_idx]:.2f}s: {e}.", rt=float(rt[peak_idx]), reason=str(e)
            )
        return emg_params

    # ==========================================================
    # 7️⃣  MÉTRICAS DO VALE
    # ==========================================================
    def find_valley(self, intensity, peak1_idx, peak2_idx):
        return peak1_idx + int(np.argmin(intensity[peak1_idx:peak2_idx]))

    def calculate_valley_metrics(self, rt, intensity, peak1_idx, peak2_idx):
        valley_idx = self.find_valley(intensity, peak1_idx, peak2_idx)
        h_valley = intensity[valley_idx]
        h1, h2 = intensity[peak1_idx], intensity[peak2_idx]
        h_menor = min(h1, h2)
        h_maior = max(h1, h2)
        valley_pct = (1.0 - h_valley / h_menor) * 100.0 if h_menor > 0 else 0.0
        height_ratio = h_menor / h_maior if h_maior > 0 else 1.0
        return valley_idx, valley_pct, height_ratio

    # ==========================================================
    # 7️⃣b CLASSIFICAÇÃO DO TIPO DE SOBREPOSIÇÃO
    # ==========================================================
    def classify_overlap(self, rt, intensity, peak1_idx, peak2_idx):
        valley_idx, valley_pct, height_ratio = self.calculate_valley_metrics(rt, intensity, peak1_idx, peak2_idx)
        m = self._m
        if valley_pct >= m.valley_pct_independent:
            return "INDEPENDENT", valley_idx, valley_pct, height_ratio, None, None
        if valley_pct >= m.valley_pct_dropline:
            return "DROP_LINE", valley_idx, valley_pct, height_ratio, None, None
        h1, h2 = intensity[peak1_idx], intensity[peak2_idx]
        if height_ratio <= m.height_ratio_rider and valley_pct <= m.valley_pct_skim_max:
            parent_idx = peak1_idx if h1 >= h2 else peak2_idx
            rider_idx = peak2_idx if h1 >= h2 else peak1_idx
            return "TANGENT_SKIM", valley_idx, valley_pct, height_ratio, parent_idx, rider_idx
        return "DECONVOLUTION", valley_idx, valley_pct, height_ratio, None, None

    # ==========================================================
    # 7️⃣c DROP-LINE
    # ==========================================================
    def integrate_dropline(self, rt, intensity, peak1_idx, peak2_idx, valley_idx, left_base, right_base):
        self.audit.info(
            "Integration",
            f"Drop-line iniciado: RT={rt[peak1_idx]:.2f}s e {rt[peak2_idx]:.2f}s | valley_idx={valley_idx}.",
            rt1=float(rt[peak1_idx]),
            rt2=float(rt[peak2_idx]),
            valley_idx=valley_idx,
        )
        results = []
        for pk_idx, l, r in [(peak1_idx, left_base, valley_idx), (peak2_idx, valley_idx, right_base)]:
            area_trap, _, _, _ = self.integrate_trapezoid_segment(rt, intensity, l, r)
            row = self.fit_emg_peak(rt, intensity, pk_idx, l, r)
            if row is None:
                y_seg = intensity[l:r]
                row = {
                    "A_param": np.nan,
                    "rt": rt[pk_idx],
                    "height": float(np.max(y_seg)) if len(y_seg) > 0 else 0.0,
                    "area": area_trap,
                    "area_emg": np.nan,
                    "sigma": np.nan,
                    "tau": np.nan,
                    "tailing_factor": np.nan,
                    "tau_bounded": np.nan,
                    "marker_rt": rt[pk_idx],
                    "marker_height": intensity[pk_idx],
                    "peak_index_apex": int(pk_idx),
                    "peak_index_start": int(l),
                    "peak_index_end": int(r),
                }
                self.audit.log_emg_fallback(
                    rt_s=float(rt[pk_idx]), area_trap=area_trap, reason="fit_emg_peak retornou None (janela drop-line)"
                )
            else:
                row["area"] = area_trap
            snr, _, noise = self.estimate_local_snr(rt, intensity, pk_idx, l, r)
            row["snr"] = snr
            row["local_noise"] = noise
            row["integration_method"] = "DROP_LINE"
            results.append(row)
        return results

    # ==========================================================
    # 7️⃣d TANGENT SKIM
    # ==========================================================
    def integrate_tangent_skim(self, rt, intensity, parent_idx, rider_idx, valley_idx, left_base_parent, right_base_rider):
        self.audit.info(
            "Integration",
            f"Tangent Skim iniciado: pai RT={rt[parent_idx]:.2f}s | rider RT={rt[rider_idx]:.2f}s.",
            parent_rt=float(rt[parent_idx]),
            rider_rt=float(rt[rider_idx]),
        )
        results = []
        rider_right = rider_idx > parent_idx
        if rider_right:
            parent_l, parent_r = left_base_parent, valley_idx
            rider_l, rider_r = valley_idx, right_base_rider
        else:
            parent_l, parent_r = valley_idx, right_base_rider
            rider_l, rider_r = left_base_parent, valley_idx
        area_parent, _, _, _ = self.integrate_trapezoid_segment(rt, intensity, parent_l, parent_r)
        row_parent = self.fit_emg_peak(rt, intensity, parent_idx, parent_l, parent_r)
        if row_parent:
            row_parent["area"] = area_parent
            snr, _, noise = self.estimate_local_snr(rt, intensity, parent_idx, parent_l, parent_r)
            row_parent["snr"] = snr
            row_parent["local_noise"] = noise
            row_parent["integration_method"] = "TANGENT_SKIM_PARENT"
            results.append(row_parent)
        x_rider = rt[rider_l : rider_r + 1]
        y_rider = intensity[rider_l : rider_r + 1]
        x0, y0 = rt[rider_l], intensity[rider_l]
        x1, y1 = rt[rider_r], intensity[rider_r]
        tangent = np.interp(x_rider, [x0, x1], [y0, y1])
        y_above = np.maximum(y_rider - tangent, 0.0)
        area_rider = float(trapezoid(y_above, x_rider))
        rider_apex = int(np.argmax(y_above))
        height_rider = float(y_above[rider_apex])
        self.audit.log_integration(
            method="TANGENT_SKIM_RIDER",
            rt=float(x_rider[rider_apex]),
            area=area_rider,
            snr=0.0,
            window=(rider_l, rider_r),
            extra=dict(tangent_y0=float(y0), tangent_y1=float(y1)),
        )
        row_rider = {
            "A_param": np.nan,
            "rt": float(x_rider[rider_apex]),
            "height": height_rider,
            "area": area_rider,
            "area_emg": np.nan,
            "sigma": np.nan,
            "tau": np.nan,
            "tailing_factor": np.nan,
            "tau_bounded": np.nan,
            "snr": np.nan,
            "local_noise": np.nan,
            "integration_method": "TANGENT_SKIM_RIDER",
            "_skim_x": x_rider,
            "_skim_tangent": tangent,
            "marker_rt": rt[rider_idx],
            "marker_height": intensity[rider_idx],
            "peak_index_apex": int(rider_l + rider_apex),
            "peak_index_start": int(rider_l),
            "peak_index_end": int(rider_r),
        }
        if len(y_above) > 4:
            noise_above = np.std(np.diff(y_above)) / np.sqrt(2)
            if noise_above > 0:
                row_rider["snr"] = height_rider / noise_above
                row_rider["local_noise"] = noise_above
        results.append(row_rider)
        return results

    # ==========================================================
    # 📐  SYSTEM SUITABILITY
    # ==========================================================
    def peak_width_at_fraction(self, rt, intensity, peak_idx, fraction: float):
        apex_height = float(intensity[peak_idx])
        threshold = fraction * apex_height
        apex_rt = float(rt[peak_idx])
        left_x = None
        for j in range(peak_idx, 0, -1):
            if intensity[j - 1] <= threshold <= intensity[j]:
                x0, y0 = rt[j - 1], intensity[j - 1]
                x1, y1 = rt[j], intensity[j]
                left_x = float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))
                break
        if left_x is None:
            left_x = float(rt[0])
        right_x = None
        for j in range(peak_idx, len(rt) - 1):
            if intensity[j] >= threshold >= intensity[j + 1]:
                x0, y0 = rt[j], intensity[j]
                x1, y1 = rt[j + 1], intensity[j + 1]
                right_x = float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))
                break
        if right_x is None:
            right_x = float(rt[-1])
        width = right_x - left_x
        front = apex_rt - left_x
        tail = right_x - apex_rt
        return width, left_x, right_x, front, tail

    def tailing_factor(self, rt, intensity, peak_idx) -> float:
        W, left_x, right_x, front, tail = self.peak_width_at_fraction(rt, intensity, peak_idx, fraction=0.05)
        return np.nan if front <= 0 else float(W / (2.0 * front))

    def asymmetry_factor(self, rt, intensity, peak_idx) -> float:
        _, left_x, right_x, front, tail = self.peak_width_at_fraction(rt, intensity, peak_idx, fraction=0.10)
        return np.nan if front <= 0 else float(tail / front)

    def theoretical_plates(self, rt, intensity, peak_idx) -> float:
        W_half, *_ = self.peak_width_at_fraction(rt, intensity, peak_idx, fraction=0.50)
        t_R = float(rt[peak_idx])
        return np.nan if W_half <= 0 or t_R <= 0 else float(5.54 * (t_R / W_half) ** 2)

    def compute_usp_metrics(self, rt: np.ndarray, intensity: np.ndarray, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calcula métricas farmacopeicas por pico e as adiciona ao DataFrame.

        Colunas adicionadas
        -------------------
        N_plates_ep       : eficiência (EP/HPLC): 5.54 × (tR / W½)²
                            usa largura a meia-altura — mais robusta para picos
                            assimétricos, exigida pela Farmacopeia Europeia.
        N_plates_usp      : eficiência (USP): 16 × (tR / Wbase)²
                            usa largura na base (5 % da altura) — fórmula clássica
                            USP; tende a ser mais afetada por tailing.
        N_plates          : alias de N_plates_ep (preserva retrocompatibilidade).
        tailing_factor_usp: T_f = W_{0.05} / (2 × d_front) — USP/JP.
        asymmetry_factor_ep: As = d_tail / d_front (a 10 % da altura) — EP.
        W_half_s          : largura do pico a 50 % da altura (segundos).
        W_base_s          : largura do pico a 5 % da altura ≈ base (segundos).
        area_pct          : participação percentual da área do pico no total.
        """
        required = {"peak_index_apex", "peak_index_start", "peak_index_end"}
        missing = required - set(df.columns)
        if missing:
            raise KeyError(f"DataFrame não contém as colunas de índice: {missing}.")
        out = df.copy()

        # ── Área total para normalização ──────────────────────────────────────
        total_area = float(out["area"].sum()) if "area" in out.columns else 0.0

        cols: dict[str, list] = {
            c: []
            for c in (
                "N_plates_ep",
                "N_plates_usp",
                "tailing_factor_usp",
                "asymmetry_factor_ep",
                "W_half_s",
                "W_base_s",
                "area_pct",
            )
        }

        for _, row in out.iterrows():
            apex = int(row["peak_index_apex"])

            if apex < 0 or apex >= len(rt):
                self.audit.warn("USPMetrics", f"peak_index_apex={apex} fora do intervalo.", apex=apex)
                for lst in cols.values():
                    lst.append(np.nan)
                continue

            # ── Larguras nas frações necessárias ─────────────────────────────
            W_half, _, _, _, _ = self.peak_width_at_fraction(rt, intensity, apex, fraction=0.50)
            W_base, _, _, Df, _ = self.peak_width_at_fraction(rt, intensity, apex, fraction=0.05)
            _, _, _, Df10, _ = self.peak_width_at_fraction(rt, intensity, apex, fraction=0.10)
            _, _, _, _, Dt = self.peak_width_at_fraction(rt, intensity, apex, fraction=0.10)

            t_R = float(rt[apex])

            # N — fórmula EP (meia-altura) ─────────────────────────────────────
            n_ep = float(5.54 * (t_R / W_half) ** 2) if W_half > 0 and t_R > 0 else np.nan
            # N — fórmula USP (base) ───────────────────────────────────────────
            n_usp = float(16.0 * (t_R / W_base) ** 2) if W_base > 0 and t_R > 0 else np.nan

            # Tailing Factor (USP) ─────────────────────────────────────────────
            tf = float(W_base / (2.0 * Df)) if Df > 0 else np.nan

            # Asymmetry Factor (EP) ────────────────────────────────────────────
            af = float(Dt / Df10) if Df10 > 0 else np.nan

            # Área % ──────────────────────────────────────────────────────────
            peak_area = float(row.get("area", np.nan))
            area_pct = float(peak_area / total_area * 100.0) if total_area > 0 and np.isfinite(peak_area) else np.nan

            cols["N_plates_ep"].append(n_ep)
            cols["N_plates_usp"].append(n_usp)
            cols["tailing_factor_usp"].append(tf)
            cols["asymmetry_factor_ep"].append(af)
            cols["W_half_s"].append(float(W_half))
            cols["W_base_s"].append(float(W_base))
            cols["area_pct"].append(area_pct)

            self.audit.info(
                "USPMetrics",
                f"RT={t_R:.2f}s: N_ep={n_ep:.0f}, N_usp={n_usp:.0f}, "
                f"TF={tf:.3f}, AF={af:.3f}, W½={W_half:.3f}s, Area%={area_pct:.2f}%.",
                rt=t_R,
                apex_idx=apex,
                N_ep=n_ep,
                N_usp=n_usp,
                tailing_factor=tf,
                asymmetry_factor=af,
                W_half_s=float(W_half),
                W_base_s=float(W_base),
                area_pct=area_pct,
            )

        for col_name, col_data in cols.items():
            out[col_name] = col_data

        # Alias retrocompatível
        out["N_plates"] = out["N_plates_ep"]
        return out

    @staticmethod
    def percent_rsd(areas: list | np.ndarray) -> dict:
        a = np.asarray(areas, dtype=float)
        a = a[np.isfinite(a)]
        if len(a) < 2:
            return {"mean": np.nan, "std": np.nan, "rsd_pct": np.nan, "n": len(a), "status": "INSUFFICIENT DATA"}
        mean = float(np.mean(a))
        std = float(np.std(a, ddof=1))
        rsd_pct = (std / mean * 100.0) if mean > 0 else np.nan
        return {
            "mean": mean,
            "std": std,
            "rsd_pct": rsd_pct,
            "n": int(len(a)),
            "status": "PASS" if (pd.notna(rsd_pct) and rsd_pct <= 2.0) else "FAIL",
        }

    # ==========================================================
    # 8️⃣  DECONVOLUÇÃO
    # ==========================================================
    def fit_overlapping_peaks(self, rt, intensity, peak_indices, left, right):
        self.audit.info(
            "Deconvolution",
            f"Deconvolução de {len(peak_indices)} picos, janela [{left}:{right}].",
            n_peaks=len(peak_indices),
            window=(left, right),
        )
        x = rt[left:right]
        y = intensity[left:right]
        num_peaks = len(peak_indices)
        area_total_trap, _, _, _ = self.integrate_trapezoid_segment(rt, intensity, left, right)
        total_area = trapezoid(y, x)
        initial_guess, initial_guess_g = [], []
        for p in peak_indices:
            A0 = total_area / num_peaks
            mu0 = rt[p]
            sigma0 = max((rt[right] - rt[left]) / (6 * num_peaks), 0.01)
            tau0 = sigma0
            initial_guess.extend([A0, mu0, sigma0, tau0])
            initial_guess_g.extend([intensity[p], mu0, sigma0])
        lower_bounds = [0, min(x), 0.001, 0.001] * num_peaks
        upper_bounds = [np.inf, max(x), np.inf, np.inf] * num_peaks
        lower_bounds_g = [0, min(x), 0.001] * num_peaks
        upper_bounds_g = [np.inf, max(x), np.inf] * num_peaks
        try:
            popt, _ = curve_fit(self.multi_emg, x, y, p0=initial_guess, bounds=(lower_bounds, upper_bounds), maxfev=10000)
        except Exception as e:
            self.audit.error("Deconvolution", f"Ajuste Multi-EMG falhou: {e}.", reason=str(e))
            return None
        popt_g = None
        try:
            popt_g, _ = curve_fit(
                self.multi_gaussian, x, y, p0=initial_guess_g, bounds=(lower_bounds_g, upper_bounds_g), maxfev=10000
            )
        except Exception as e:
            self.audit.warn("Deconvolution", f"Ajuste Multi-Gaussiano falhou: {e}.", reason=str(e))
        emg_areas = []
        for i in range(num_peaks):
            A, mu, sigma, tau = popt[i * 4], popt[i * 4 + 1], abs(popt[i * 4 + 2]), abs(popt[i * 4 + 3])
            emg_areas.append(float(trapezoid(self.emg(x, A, mu, sigma, tau), x)))
        total_emg_area = sum(emg_areas)
        self.audit.log_deconv_audit(area_total_trap, total_emg_area)
        results = []
        for i in range(num_peaks):
            A, mu, sigma, tau = popt[i * 4], popt[i * 4 + 1], abs(popt[i * 4 + 2]), abs(popt[i * 4 + 3])
            y_comp = self.emg(x, A, mu, sigma, tau)
            frac_emg = emg_areas[i] / total_emg_area if total_emg_area > 0 else 1.0 / num_peaks
            area_comp_trap = frac_emg * area_total_trap
            self.audit.info(
                "Deconvolution",
                f"Componente {i+1}: RT={mu:.2f}s, fração EMG={frac_emg:.3f}, área={area_comp_trap:.0f}.",
                component=i + 1,
                rt=float(mu),
                frac_emg=frac_emg,
                area=area_comp_trap,
            )
            row = {
                "A_param": A,
                "rt": mu,
                "height": float(np.max(y_comp)),
                "area": area_comp_trap,
                "area_emg": emg_areas[i],
                "area_frac_emg": frac_emg,
                "sigma": sigma,
                "tau": tau,
                "tailing_factor": tau / sigma if sigma > 0 else np.nan,
                "integration_method": "DECONVOLUTION",
                "marker_rt": rt[peak_indices[i]],
                "marker_height": intensity[peak_indices[i]],
                "peak_index_apex": int(peak_indices[i]),
                "peak_index_start": int(left),
                "peak_index_end": int(right),
            }
            if popt_g is not None:
                row.update({"gauss_A": popt_g[i * 3], "gauss_mu": popt_g[i * 3 + 1], "gauss_sigma": abs(popt_g[i * 3 + 2])})
            results.append(row)
        return results

    # ==========================================================
    # 9️⃣  REMOÇÃO DE PICO DO SOLVENTE
    # ==========================================================
    def remove_solvent_peak(self, df, protected_rt: Optional[float] = None):
        """
        Remove picos de solvente do DataFrame de integração.

        Parâmetros
        ----------
        df : pd.DataFrame
            DataFrame retornado por ``integrate()`` antes da limpeza.
        protected_rt : float | None
            RT (em segundos) do Padrão Interno detectado ANTES desta remoção.
            Se fornecido, qualquer pico dentro da janela
            [protected_rt ± is_search_window_s] é preservado mesmo que
            satisfaça os critérios de remoção de solvente.

            Motivação: um IS de grande área (ex.: solvente deuterado ou padrão
            interno majoritário) seria erroneamente descartado pelo filtro
            ``area > factor × mediana`` por ter área muito maior que os analitos.
            Protegê-lo aqui garante que ``align_runs()`` sempre o encontre.
        """
        if df.empty:
            return df
        rt_min_exclude = self._m.solvent_rt_cutoff_s
        area_factor = self._m.solvent_area_factor
        median_area = float(np.median(df["area"]))

        keep_mask = (df["rt"] > rt_min_exclude) & (df["area"] < area_factor * median_area)

        # ── ✅ FIX 1: Protege o IS da remoção por área ───────────────────────
        # Problema original: IS com RT > solvent_rt_cutoff_s mas área >> mediana
        # era removido pelo filtro de área. Ao receber o RT do IS pré-identificado,
        # adicionamos uma cláusula OR que força a preservação do pico independente
        # de sua área.
        if protected_rt is not None:
            is_window = self._m.is_search_window_s
            is_mask = (df["rt"] - protected_rt).abs() <= is_window
            if is_mask.any():
                self.audit.info(
                    "QC",
                    f"IS protegido da remoção de solvente: RT≈{protected_rt:.2f}s "
                    f"(janela ±{is_window:.1f}s). "
                    f"{int(is_mask.sum())} pico(s) preservado(s) independente da área.",
                    is_rt=float(protected_rt),
                    is_window_s=float(is_window),
                    n_protected=int(is_mask.sum()),
                )
            keep_mask = keep_mask | is_mask

        filtered = df[keep_mask]
        self.audit.log_solvent_removal(
            n_before=len(df),
            n_after=len(filtered),
            median_area=median_area,
            rt_cutoff=rt_min_exclude,
            factor=area_factor,
        )
        return filtered.reset_index(drop=True)

    # ==========================================================
    # 🔟  PIPELINE DE INTEGRAÇÃO COMPLETO
    # ==========================================================
    def integrate(self, rt, intensity):
        """
        Executa o pipeline completo de integração.
        Todos os parâmetros vêm de self.method.
        """
        snr_threshold = self._m.snr_threshold
        rs_threshold = self._m.rs_deconv_threshold
        self.audit.info("Pipeline", "Integração iniciada.", snr_threshold=snr_threshold)
        peaks, left_ips, right_ips, snr_values, left_base, right_base = self.find_peaks(rt, intensity)
        results = []
        i = 0
        while i < len(peaks):
            if i < len(peaks) - 1:
                t1 = rt[peaks[i]]
                t2 = rt[peaks[i + 1]]
                dt = rt[1] - rt[0]
                w1 = (right_ips[i] - left_ips[i]) * dt
                w2 = (right_ips[i + 1] - left_ips[i + 1]) * dt
                Rs = 2 * (t2 - t1) / (w1 + w2) if (w1 + w2) > 0 else 99.0
                if Rs < rs_threshold:
                    left_global = int(left_base[i])
                    right_global = int(right_base[i + 1])
                    method, valley_idx, valley_pct, height_ratio, parent_idx, rider_idx = self.classify_overlap(
                        rt, intensity, peaks[i], peaks[i + 1]
                    )
                    self.audit.log_overlap_decision(
                        rt1=float(t1),
                        rt2=float(t2),
                        Rs=float(Rs),
                        valley_pct=float(valley_pct),
                        height_ratio=float(height_ratio),
                        method=method,
                    )
                    if method == "INDEPENDENT":
                        overlap = []
                        for j in range(2):
                            pk_idx = peaks[i + j]
                            row = self.fit_emg_peak(rt, intensity, pk_idx, int(left_base[i + j]), int(right_base[i + j]))
                            if row:
                                snr, _, noise = self.estimate_local_snr(
                                    rt, intensity, pk_idx, int(left_base[i + j]), int(right_base[i + j])
                                )
                                row["snr"] = snr
                                row["local_noise"] = noise
                                row["integration_method"] = "EMG"
                                overlap.append(row)
                    elif method == "DROP_LINE":
                        overlap = self.integrate_dropline(
                            rt, intensity, peaks[i], peaks[i + 1], valley_idx, left_global, right_global
                        )
                    elif method == "TANGENT_SKIM":
                        overlap = self.integrate_tangent_skim(
                            rt, intensity, parent_idx, rider_idx, valley_idx, left_global, right_global
                        )
                    else:
                        overlap = self.fit_overlapping_peaks(rt, intensity, [peaks[i], peaks[i + 1]], left_global, right_global)
                        if overlap:
                            for j, ov_row in enumerate(overlap):
                                pk_idx = peaks[i + j]
                                snr, _, noise = self.estimate_local_snr(
                                    rt, intensity, pk_idx, int(left_base[i + j]), int(right_base[i + j])
                                )
                                ov_row["snr"] = snr
                                ov_row["local_noise"] = noise
                    if overlap:
                        for ov_row in overlap:
                            ov_row["valley_pct"] = valley_pct
                            ov_row["height_ratio"] = height_ratio
                            ov_row["Rs"] = Rs
                        results.extend(overlap)
                        self.audit.info(
                            "Integration",
                            f"{len(overlap)} pico(s) adicionado(s) via {method} (RT≈{t1:.1f}–{t2:.1f}s).",
                            method=method,
                            n_peaks=len(overlap),
                            rt1=float(t1),
                            rt2=float(t2),
                        )
                        i += 2
                        continue
            peak_result = self.fit_emg_peak(rt, intensity, peaks[i], int(left_base[i]), int(right_base[i]))
            if peak_result:
                snr, _, noise = self.estimate_local_snr(rt, intensity, peaks[i], int(left_base[i]), int(right_base[i]))
                peak_result["snr"] = snr
                peak_result["local_noise"] = noise
                peak_result["integration_method"] = "EMG"
                peak_result["valley_pct"] = np.nan
                peak_result["height_ratio"] = np.nan
                peak_result["Rs"] = np.nan
                results.append(peak_result)
                self.audit.info(
                    "Integration",
                    f"Pico isolado integrado: RT={peak_result['rt']:.2f}s, área={peak_result['area']:.0f}, SNR={snr:.1f}.",
                    rt=peak_result["rt"],
                    area=peak_result["area"],
                    snr=snr,
                )
            i += 1

        self._skim_traces = []
        clean_results = []
        for row in results:
            skim_x = row.pop("_skim_x", None)
            skim_tan = row.pop("_skim_tangent", None)
            self._skim_traces.append((skim_x, skim_tan) if skim_x is not None else None)
            clean_results.append(row)

        df = pd.DataFrame(clean_results)

        # ── ✅ FIX 2: Detecta o IS ANTES de remove_solvent_peak ──────────────
        # Problema original: integrate() chamava remove_solvent_peak() sem passar
        # o RT do IS. Se o IS tinha área >> mediana (ex.: solvente deuterado com
        # área 50× maior que os analitos), o filtro `area > factor × mediana`
        # o descartava silenciosamente. align_runs() recebia então um DataFrame
        # sem o IS → find_internal_standard() falhava → ValueError → corrida
        # ignorada → enriched vazio → ValueError fatal.
        #
        # Correção: identificamos o IS aqui, passamos seu RT para
        # remove_solvent_peak() que o protege explicitamente via OR-mask.
        protected_is_rt: Optional[float] = None
        if not df.empty and self._m.is_rt_seconds is not None:
            try:
                is_row = self.find_internal_standard(df)
                protected_is_rt = float(is_row["rt"])
                self.audit.info(
                    "Pipeline",
                    f"IS pré-identificado antes da remoção de solvente: "
                    f"RT={protected_is_rt:.2f}s — será protegido do filtro de área.",
                    is_rt=protected_is_rt,
                )
            except ValueError as exc:
                self.audit.warn(
                    "Pipeline",
                    f"IS não localizado antes da remoção de solvente ({exc}). " "Proteção do IS desativada nesta corrida.",
                    reason=str(exc),
                )

        df = self.remove_solvent_peak(df, protected_rt=protected_is_rt)

        self.audit.info("Pipeline", f"Integração concluída: {len(df)} pico(s) no relatório final.", n_peaks_final=len(df))
        return df

    # ==========================================================
    # 🔗  ALINHAMENTO MULTI-CORRIDA — RRT + BINNING
    # ==========================================================

    def find_internal_standard(self, df: pd.DataFrame) -> pd.Series:
        if df.empty:
            raise ValueError("DataFrame vazio — não é possível localizar o IS.")
        m = self._m
        if m.is_rt_seconds is not None:
            lo = m.is_rt_seconds - m.is_search_window_s
            hi = m.is_rt_seconds + m.is_search_window_s
            window = df[(df["rt"] >= lo) & (df["rt"] <= hi)]
            if window.empty:
                raise ValueError(
                    f"Nenhum pico encontrado na janela IS [{lo:.1f}–{hi:.1f}]s. "
                    f"Verifique is_rt_seconds ({m.is_rt_seconds}s) e is_search_window_s ({m.is_search_window_s}s)."
                )
            is_row = window.loc[window["area"].idxmax()]
            strategy = f"maior área em [{lo:.1f}–{hi:.1f}]s"
        else:
            is_row = df.loc[df["area"].idxmax()]
            strategy = "maior área global (fallback automático)"
        self.audit.info(
            "RRTAlignment",
            f"IS localizado: RT={is_row['rt']:.3f}s, área={is_row['area']:.0f} ({strategy}).",
            is_rt=float(is_row["rt"]),
            is_area=float(is_row["area"]),
            strategy=strategy,
        )
        return is_row

    def compute_rrt(self, df: pd.DataFrame, is_rt: float) -> pd.DataFrame:
        if is_rt <= 0:
            raise ValueError(f"is_rt deve ser positivo; recebido: {is_rt}.")
        out = df.copy()
        out["rrt"] = out["rt"] / is_rt
        out["is_rt_used"] = float(is_rt)
        self.audit.info(
            "RRTAlignment",
            f"RRT calculado para {len(out)} picos usando IS RT={is_rt:.3f}s.",
            is_rt=float(is_rt),
            n_peaks=len(out),
            rrt_min=float(out["rrt"].min()),
            rrt_max=float(out["rrt"].max()),
        )
        return out

    def align_runs(self, runs: list[tuple[str, pd.DataFrame]], auto_is: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
        if not runs:
            raise ValueError("Lista de corridas vazia.")
        m = self._m
        tol = m.rrt_bin_tolerance
        self.audit.info(
            "RRTAlignment",
            f"align_runs iniciado: {len(runs)} corridas, tolerância RRT={tol}.",
            n_runs=len(runs),
            rrt_bin_tolerance=tol,
            auto_is=auto_is,
        )

        enriched: list[pd.DataFrame] = []
        prev_is_rt: Optional[float] = None  # rastreia IS da corrida anterior para detectar deriva

        for run_id, df in runs:
            if auto_is:
                try:
                    is_row = self.find_internal_standard(df)
                    is_rt = float(is_row["rt"])
                except ValueError as exc:
                    self.audit.error(
                        "RRTAlignment",
                        f"Corrida '{run_id}': IS não encontrado — {exc}. Corrida ignorada.",
                        run_id=run_id,
                        error=str(exc),
                    )
                    continue

                # ── ✅ FIX 3: Detecta deriva significativa do IS entre corridas ──
                # Emite WARN quando o RT do IS muda mais que is_search_window_s
                # entre corridas consecutivas. Isso pode indicar: (a) o IS foi
                # substituído por um composto diferente, (b) deriva cromatográfica
                # grave, ou (c) erro de configuração no is_rt_seconds.
                if prev_is_rt is not None:
                    delta_is = abs(is_rt - prev_is_rt)
                    if delta_is > m.is_search_window_s:
                        self.audit.warn(
                            "RRTAlignment",
                            f"Corrida '{run_id}': RT do IS deslocou {delta_is:.2f}s em relação à corrida anterior "
                            f"(anterior={prev_is_rt:.2f}s, atual={is_rt:.2f}s). "
                            f"Deslocamento > is_search_window_s ({m.is_search_window_s:.1f}s). "
                            "Verifique condições cromatográficas ou identidade do IS.",
                            original_value=round(prev_is_rt, 3),
                            new_value=round(is_rt, 3),
                            run_id=run_id,
                            delta_is_s=round(delta_is, 3),
                            is_search_window_s=m.is_search_window_s,
                        )
                prev_is_rt = is_rt

                rrt_df = self.compute_rrt(df, is_rt)
            else:
                if "rrt" not in df.columns:
                    self.audit.error(
                        "RRTAlignment",
                        f"Corrida '{run_id}' não contém coluna 'rrt' e auto_is=False. Corrida ignorada.",
                        run_id=run_id,
                    )
                    continue
                rrt_df = df.copy()
            rrt_df["run_id"] = run_id
            enriched.append(rrt_df)

        if not enriched:
            raise ValueError("Nenhuma corrida pôde ser processada. Verifique os parâmetros do IS e os DataFrames de entrada.")

        all_peaks = pd.concat(enriched, ignore_index=True).sort_values("rrt").reset_index(drop=True)
        bin_ids: list[int] = []
        current_bin = -1
        current_centroid = -np.inf
        current_count = 0
        for rrt_val in all_peaks["rrt"]:
            if abs(rrt_val - current_centroid) <= tol:
                current_centroid = (current_centroid * current_count + rrt_val) / (current_count + 1)
                current_count += 1
            else:
                current_bin += 1
                current_centroid = float(rrt_val)
                current_count = 1
            bin_ids.append(current_bin)
        all_peaks["bin_id"] = bin_ids
        all_peaks["bin_rrt_centroid"] = [all_peaks.loc[all_peaks["bin_id"] == b, "rrt"].mean() for b in bin_ids]
        n_bins = all_peaks["bin_id"].nunique()
        self.audit.info(
            "RRTAlignment",
            f"Binning concluído: {len(all_peaks)} picos em {n_bins} bins (tol={tol}).",
            n_peaks_total=len(all_peaks),
            n_bins=n_bins,
        )

        cols_long = [
            "bin_id",
            "bin_rrt_centroid",
            "run_id",
            "rt",
            "rrt",
            "area",
            "integration_method",
            "snr",
            "peak_index_apex",
            "peak_index_start",
            "peak_index_end",
        ]
        df_long = all_peaks[[c for c in cols_long if c in all_peaks.columns]].copy()

        stat_rows = []
        n_runs_total = len(runs)
        for bin_id, group in all_peaks.groupby("bin_id"):
            centroid = float(group["bin_rrt_centroid"].iloc[0])
            n_present = group["run_id"].nunique()
            areas = group["area"].dropna().values
            rrt_mean = float(group["rrt"].mean())
            rrt_std = float(group["rrt"].std(ddof=1)) if len(group) > 1 else 0.0
            rrt_cv = (rrt_std / rrt_mean * 100.0) if rrt_mean > 0 else np.nan
            area_mean = float(np.mean(areas)) if len(areas) > 0 else np.nan
            area_std = float(np.std(areas, ddof=1)) if len(areas) > 1 else 0.0
            area_cv = (area_std / area_mean * 100.0) if (area_mean and area_mean > 0) else np.nan
            area_status = "PASS" if (pd.notna(area_cv) and area_cv <= 2.0) else "FAIL" if pd.notna(area_cv) else "N/A"
            self.audit.info(
                "RRTAlignment",
                f"Bin {bin_id}: centróide RRT={centroid:.4f}, {n_present}/{n_runs_total} corridas, CV(área)={area_cv:.1f}%.",
                bin_id=int(bin_id),
                rrt_centroid=centroid,
                n_runs_present=n_present,
                area_cv_pct=area_cv,
                area_status=area_status,
            )
            stat_rows.append(
                {
                    "bin_id": int(bin_id),
                    "bin_rrt_centroid": centroid,
                    "n_runs_present": int(n_present),
                    "n_runs_total": int(n_runs_total),
                    "rrt_mean": rrt_mean,
                    "rrt_cv_pct": rrt_cv,
                    "area_mean": area_mean,
                    "area_cv_pct": area_cv,
                    "area_status": area_status,
                    "rt_mean_s": float(group["rt"].mean()),
                    "rt_std_s": float(group["rt"].std(ddof=1)) if len(group) > 1 else 0.0,
                }
            )
        df_stats = pd.DataFrame(stat_rows)
        self.audit.info(
            "RRTAlignment",
            f"align_runs concluído. Bins totais: {n_bins} | completos: {int((df_stats['n_runs_present'] == n_runs_total).sum())}.",
            n_bins=n_bins,
            n_bins_complete=int((df_stats["n_runs_present"] == n_runs_total).sum()),
        )
        return df_long, df_stats

    # ==========================================================
    # 🔄  PROCESSAMENTO EM LOTE
    # ==========================================================
    def process_batch(self, cdf_files: list[str | Path], *, compute_usp: bool = True, align: bool = False):
        n_total = len(cdf_files)
        self.audit.info(
            "Batch", f"Lote iniciado: {n_total} arquivo(s) CDF.", n_files=n_total, compute_usp=compute_usp, align=align
        )
        results: list[RunResult] = []
        for file_idx, cdf_path in enumerate(cdf_files, start=1):
            cdf_path = Path(cdf_path)
            run_id = cdf_path.stem
            self.audit.info("Batch", f"[{file_idx}/{n_total}] Iniciando: {cdf_path.name}", run_id=run_id, file_idx=file_idx)
            rt = intensity = corrected = baseline = df = None
            try:
                rt, intensity = self.read_cdf(str(cdf_path))
                corrected, baseline = self.remove_baseline(rt, intensity)
                df = self.integrate(rt, corrected)
                if df.empty:
                    raise IntegrationError("integrate() retornou DataFrame vazio.", context={"run_id": run_id})
                if compute_usp:
                    df = self.compute_usp_metrics(rt, corrected, df)
                self.audit.info(
                    "Batch", f"[{file_idx}/{n_total}] OK: {len(df)} pico(s) — {cdf_path.name}", run_id=run_id, n_peaks=len(df)
                )
                results.append(
                    RunResult(
                        run_id=run_id, status="OK", cdf_path=str(cdf_path), results_df=df, audit_events=self.audit.to_dict_list()
                    )
                )
            except CDFReadError as exc:
                self._handle_batch_error(results, run_id, cdf_path, exc, file_idx, n_total, "Arquivo CDF ilegível.")
            except BaselineError as exc:
                self._handle_batch_error(results, run_id, cdf_path, exc, file_idx, n_total, "Falha na subtração de baseline.")
            except PeakDetectionError as exc:
                self._handle_batch_error(results, run_id, cdf_path, exc, file_idx, n_total, "Nenhum pico detectado.")
            except IntegrationError as exc:
                self._handle_batch_error(results, run_id, cdf_path, exc, file_idx, n_total, "Falha na integração.")
            except GCAnalyzerError as exc:
                self._handle_batch_error(results, run_id, cdf_path, exc, file_idx, n_total, "Erro interno do GCAnalyzer.")
            except MemoryError as exc:
                self._handle_batch_error(results, run_id, cdf_path, exc, file_idx, n_total, "Memória insuficiente.")
            except Exception as exc:
                self._handle_batch_error(results, run_id, cdf_path, exc, file_idx, n_total, "Erro inesperado.")
            finally:
                del rt, intensity, corrected, baseline, df
                _gc.collect()
        n_ok = sum(1 for r in results if r.ok)
        n_fail = len(results) - n_ok
        self.audit.info(
            "Batch",
            f"Lote concluído: {n_ok} OK | {n_fail} falha(s) | {n_total} total.",
            n_ok=n_ok,
            n_fail=n_fail,
            n_total=n_total,
        )
        df_stats: pd.DataFrame | None = None
        if align:
            ok_runs = [(r.run_id, r.results_df) for r in results if r.ok]
            if len(ok_runs) >= 2:
                try:
                    _, df_stats = self.align_runs(ok_runs)
                    self.audit.info("Batch", f"Alinhamento concluído sobre {len(ok_runs)} corridas OK.", n_aligned=len(ok_runs))
                except AlignmentError as exc:
                    self.audit.error("Batch", f"Alinhamento falhou após lote: {exc}.", error=str(exc))
            else:
                self.audit.warn(
                    "Batch",
                    f"Alinhamento solicitado mas apenas {len(ok_runs)} corrida(s) OK (mínimo: 2).",
                    n_ok_runs=len(ok_runs),
                )
        return results, df_stats

    def _handle_batch_error(self, results, run_id, cdf_path, exc, file_idx, n_total, summary):
        tb_str = traceback.format_exc()
        exc_type = type(exc).__name__
        self.audit.error(
            "Batch",
            f"[{file_idx}/{n_total}] FAILED — {run_id}: {summary} ({exc_type}: {exc})",
            run_id=run_id,
            cdf_path=str(cdf_path),
            error_type=exc_type,
            error_message=str(exc),
            traceback=tb_str,
        )
        results.append(
            RunResult(
                run_id=run_id,
                status="FAILED",
                cdf_path=str(cdf_path),
                results_df=None,
                audit_events=self.audit.to_dict_list(),
                error_type=exc_type,
                error_message=str(exc),
                error_traceback=tb_str,
            )
        )

    # ==========================================================
    # 📊  MÉTRICAS ESTENDIDAS POR PICO
    # ==========================================================
    def compute_extended_metrics(
        self,
        rt: np.ndarray,
        intensity: np.ndarray,
        df: pd.DataFrame,
        dead_time_s: Optional[float] = None,
    ) -> pd.DataFrame:
        """
        Adiciona métricas cromatográficas avançadas ao DataFrame de picos.

        Deve ser chamado **após** ``compute_usp_metrics`` (que fornece
        ``W_half_s``, ``W_base_s`` e ``tailing_factor_usp``).

        Colunas adicionadas
        -------------------
        Rs_usp
            Resolução (USP): ``Rs = 2(tR₂−tR₁) / (Wb₁+Wb₂)``
            Usa largura na base (5 % da altura ≈ 4σ).
            NaN para o primeiro pico (sem predecessor).
        Rs_ep
            Resolução (EP): ``Rs = 1.18(tR₂−tR₁) / (W½₁+W½₂)``
            Usa largura a meia-altura — menos sensível a tailing.
            NaN para o primeiro pico.
        k_prime
            Fator de retenção/capacidade: ``k' = (tR − t₀) / t₀``
            Requer ``dead_time_s > 0`` ou ``method.dead_time_s > 0``.
            NaN se t₀ não disponível.
        alpha
            Seletividade: ``α = k'ᵢ / k'ᵢ₋₁``
            Razão dos fatores de capacidade de picos consecutivos.
            NaN para o primeiro pico ou se k' não disponível.
        rsd_pct_area
            RSD% da área deste pico dentro da corrida atual (calculado sobre
            todos os picos com o mesmo bin_id, se disponível; caso contrário,
            calculado sobre todas as áreas da corrida — útil como indicador
            de homogeneidade interna do cromatograma).
        shape_quality_score (SQS)
            Score de qualidade de forma 0–1.
            Baseia-se na razão EMG τ/σ (``tailing_factor`` coluna EMG):
              - SQS = 1.0  → pico Gaussiano puro (τ/σ ≤ 1)
              - SQS → 0    → tailing severo ou adsorção na coluna
            Fórmula: ``SQS = exp(−max(τ/σ − 1, 0))``
            Fallback para ``tailing_factor_usp`` quando o EMG não convergiu.
        CQI
            Chromatographic Quality Index — score composto 0–1 por pico.
            Média geométrica ponderada de quatro sub-scores normalizados,
            controlados pelos parâmetros ``cqi_weight_*`` e ``cqi_*_ref``
            do ``ProcessingMethod``:

              N_score   = min(N_plates_ep / cqi_n_ref, 1.0)
              Rs_score  = min(Rs_usp / cqi_rs_ref, 1.0)   [0.5 se Rs NaN]
              TF_score  = max(0, 1 − |TF_usp − 1| / 2)   [penaliza afastamento de 1]
              SNR_score = min(snr / cqi_snr_ref, 1.0)

            ``CQI = (N_score^w_N × Rs_score^w_Rs × TF_score^w_TF × SNR_score^w_SNR)
                    ^ (1/(w_N+w_Rs+w_TF+w_SNR))``

        Parameters
        ----------
        rt, intensity : np.ndarray
            Mesmos arrays usados em ``integrate()``.
        df : pd.DataFrame
            DataFrame retornado por ``compute_usp_metrics()``.
        dead_time_s : float | None
            Sobrescreve ``method.dead_time_s`` se fornecido.
            0 ou None → k' e α ficam NaN.

        Returns
        -------
        pd.DataFrame
            Cópia do DataFrame com as colunas acima adicionadas.
        """
        m = self._m
        t0 = dead_time_s if dead_time_s is not None else m.dead_time_s
        if t0 <= 0:
            t0 = None  # desativa k' / α

        out = df.reset_index(drop=True).copy()

        # ── Pré-calcular W_half_s e W_base_s se ausentes ─────────────────────
        if "W_half_s" not in out.columns or "W_base_s" not in out.columns:
            self.audit.warn(
                "ExtendedMetrics",
                "W_half_s/W_base_s ausentes — chamando compute_usp_metrics internamente.",
            )
            out = self.compute_usp_metrics(rt, intensity, out)

        rs_usp_col, rs_ep_col, kp_col, alpha_col, sqs_col, cqi_col, rsd_col = ([], [], [], [], [], [], [])

        n = len(out)
        for i, row in out.iterrows():
            tR = float(row["rt"]) if "rt" in row.index else np.nan

            # ── Rs USP e Rs EP (vs pico anterior) ────────────────────────────
            # Convenção: Rs do pico i é calculado em relação ao pico i-1.
            # O primeiro pico recebe NaN (sem predecessor).
            if i == 0:
                rs_usp_col.append(np.nan)
                rs_ep_col.append(np.nan)
            else:
                prev = out.iloc[i - 1]
                tR_prev = float(prev["rt"]) if "rt" in prev.index else np.nan
                Wb = float(row.get("W_base_s", np.nan))
                Wb_prev = float(prev.get("W_base_s", np.nan))
                Wh = float(row.get("W_half_s", np.nan))
                Wh_prev = float(prev.get("W_half_s", np.nan))

                delta = tR - tR_prev
                rs_usp = (
                    (2.0 * delta / (Wb_prev + Wb))
                    if (np.isfinite(Wb) and np.isfinite(Wb_prev) and (Wb_prev + Wb) > 0)
                    else np.nan
                )
                rs_ep = (
                    (1.18 * delta / (Wh_prev + Wh))
                    if (np.isfinite(Wh) and np.isfinite(Wh_prev) and (Wh_prev + Wh) > 0)
                    else np.nan
                )
                rs_usp_col.append(float(rs_usp))
                rs_ep_col.append(float(rs_ep))

            # ── k' e α ───────────────────────────────────────────────────────
            if t0 is not None and np.isfinite(tR) and tR > t0:
                kp = (tR - t0) / t0
            else:
                kp = np.nan
            kp_col.append(kp)

            if i > 0 and np.isfinite(kp) and np.isfinite(kp_col[i - 1]) and kp_col[i - 1] > 0:
                alpha_col.append(kp / kp_col[i - 1])
            else:
                alpha_col.append(np.nan)

            # ── Shape Quality Score (SQS) ─────────────────────────────────────
            # Usa a razão τ/σ do ajuste EMG quando disponível.
            # A coluna "tailing_factor" (EMG, sem sufixo) armazena τ/σ.
            # "tailing_factor_usp" é a versão geométrica (W/2F a 5%).
            emg_tf = row.get("tailing_factor", np.nan)  # τ/σ do EMG
            if pd.notna(emg_tf) and np.isfinite(emg_tf):
                sqs = float(np.exp(-max(float(emg_tf) - 1.0, 0.0)))
            else:
                # fallback via TF geométrico USP
                usp_tf = row.get("tailing_factor_usp", np.nan)
                if pd.notna(usp_tf) and np.isfinite(usp_tf):
                    sqs = float(max(0.0, min(1.0, 2.0 - float(usp_tf))))
                else:
                    sqs = np.nan
            sqs_col.append(sqs)

            # ── RSD% de área (dentro da corrida) ─────────────────────────────
            # Calculado sobre todas as áreas do DataFrame — é um indicador da
            # homogeneidade relativa dos picos nesta corrida.
            rsd_col.append(np.nan)  # preenchido em bloco abaixo

            # ── CQI ──────────────────────────────────────────────────────────
            # Sub-scores normalizados
            N_val = float(row.get("N_plates_ep", np.nan))
            snr_val = float(row.get("snr", np.nan))
            tf_val = float(row.get("tailing_factor_usp", np.nan))
            rs_val = rs_usp_col[-1]

            n_score = min(N_val / m.cqi_n_ref, 1.0) if np.isfinite(N_val) else 0.0
            snr_score = min(snr_val / m.cqi_snr_ref, 1.0) if np.isfinite(snr_val) else 0.0
            tf_score = max(0.0, 1.0 - abs(tf_val - 1.0) / 2.0) if np.isfinite(tf_val) else 0.0
            # Rs: 0.5 para primeiro pico (sem referência), capped a 1.0
            rs_score = min(float(rs_val) / m.cqi_rs_ref, 1.0) if np.isfinite(rs_val) else 0.5

            w = (m.cqi_weight_n, m.cqi_weight_rs, m.cqi_weight_tf, m.cqi_weight_snr)
            scores = (n_score, rs_score, tf_score, snr_score)
            w_total = sum(w)
            if w_total > 0 and all(np.isfinite(s) for s in scores):
                # média geométrica ponderada
                log_sum = sum(wi * np.log(max(s, 1e-9)) for wi, s in zip(w, scores))
                cqi_val = float(np.exp(log_sum / w_total))
            else:
                cqi_val = np.nan
            cqi_col.append(cqi_val)

        # ── RSD% de área calculado em bloco ──────────────────────────────────
        areas = out["area"].values.astype(float)
        finite_areas = areas[np.isfinite(areas)]
        if len(finite_areas) >= 2:
            rsd_val = float(np.std(finite_areas, ddof=1) / np.mean(finite_areas) * 100.0)
        else:
            rsd_val = np.nan
        rsd_col = [rsd_val] * len(out)

        out["Rs_usp"] = rs_usp_col
        out["Rs_ep"] = rs_ep_col
        out["k_prime"] = kp_col
        out["alpha"] = alpha_col
        out["shape_quality_score"] = sqs_col
        out["rsd_pct_area"] = rsd_col
        out["CQI"] = cqi_col

        self.audit.info(
            "ExtendedMetrics",
            f"Métricas estendidas calculadas para {len(out)} pico(s). "
            f"dead_time_s={'OFF' if t0 is None else f'{t0:.2f}s'}, "
            f"CQI_médio={np.nanmean(cqi_col):.3f}.",
            n_peaks=len(out),
            dead_time_s=t0,
            cqi_mean=float(np.nanmean(cqi_col)) if cqi_col else np.nan,
            rsd_pct_area=rsd_val,
        )
        return out

    # ==========================================================
    # 🌍  MÉTRICAS GLOBAIS DA CORRIDA
    # ==========================================================
    def compute_global_metrics(
        self,
        rt: np.ndarray,
        raw: np.ndarray,
        corrected: np.ndarray,
        baseline: np.ndarray,
        df: pd.DataFrame,
    ) -> dict:
        """
        Calcula métricas que caracterizam a corrida **como um todo**,
        independentemente de picos individuais.

        Retorna um dict com as seguintes chaves
        -----------------------------------------
        baseline_drift
            Diferença absoluta entre o valor médio do baseline no último
            e no primeiro décimo da corrida (unidades de intensidade).
            Valores altos indicam coluna mal condicionada, vazamento de
            temperatura, ou rampa de solvente não estabilizada.
        baseline_drift_pct
            Drift normalizado pelo valor médio absoluto do baseline × 100.
        baseline_noise_sigma
            Desvio padrão do baseline (MAD × 1.4826) — mede a estabilidade
            da linha de base, não apenas a deriva direcional.
        global_snr
            SNR global da corrida: pico máximo do sinal corrigido dividido
            pelo ruído estimado do baseline.
        total_integrated_area
            Soma de todas as áreas de picos no DataFrame (unidades × segundos).
        area_coverage_pct
            Percentual do sinal bruto total que está "dentro de picos"
            (razão entre a área integrada e a área total sob o sinal corrigido).
        n_peaks
            Número de picos no DataFrame final (após filtros de solvente, etc.).
        mean_cqi
            Média do CQI de todos os picos (NaN se a coluna não existir).
        overall_quality_score (OQS)
            Score composto 0–1 para toda a corrida:
            ``OQS = CQI_médio × baseline_score × snr_score``
            onde:
              - baseline_score = max(0, 1 − drift_pct/50) [drift de 50% = score 0]
              - snr_score      = min(global_snr / 10, 1.0)
        run_duration_s
            Duração total da corrida em segundos (último RT − primeiro RT).
        sampling_rate_hz
            Taxa de amostragem estimada (pontos / segundo).

        Parameters
        ----------
        rt, raw, corrected, baseline : np.ndarray
            Arrays retornados por ``read_cdf`` e ``remove_baseline``.
        df : pd.DataFrame
            DataFrame de picos (saída de ``integrate`` ou ``compute_extended_metrics``).
        """
        # ── Baseline: derive e ruído ──────────────────────────────────────────
        n_pts = len(baseline)
        tenth = max(1, n_pts // 10)
        bl_start = float(np.mean(baseline[:tenth]))
        bl_end = float(np.mean(baseline[-tenth:]))
        drift_abs = abs(bl_end - bl_start)
        bl_mean_abs = float(np.mean(np.abs(baseline)))
        drift_pct = (drift_abs / bl_mean_abs * 100.0) if bl_mean_abs > 0 else 0.0

        bl_median = float(np.median(baseline))
        bl_mad = float(np.median(np.abs(baseline - bl_median)))
        bl_noise = bl_mad * 1.4826

        # ── SNR global ───────────────────────────────────────────────────────
        peak_max = float(np.max(corrected))
        global_snr = (peak_max / bl_noise) if bl_noise > 0 else np.nan

        # ── Área total integrada ──────────────────────────────────────────────
        total_integrated = float(df["area"].sum()) if not df.empty and "area" in df.columns else 0.0

        # ── Coverage: % do sinal que está em picos ────────────────────────────
        total_signal_area = float(trapezoid(np.maximum(corrected, 0.0), rt))
        area_coverage = (total_integrated / total_signal_area * 100.0) if total_signal_area > 0 else 0.0

        # ── CQI médio ─────────────────────────────────────────────────────────
        mean_cqi = float(df["CQI"].mean()) if "CQI" in df.columns and not df.empty else np.nan

        # ── Overall Quality Score ─────────────────────────────────────────────
        baseline_score = max(0.0, 1.0 - drift_pct / 50.0)
        snr_score_g = min(global_snr / 10.0, 1.0) if np.isfinite(global_snr) else 0.0
        cqi_score = mean_cqi if np.isfinite(mean_cqi) else 0.5
        oqs = float(cqi_score * baseline_score * snr_score_g)

        # ── Meta da corrida ───────────────────────────────────────────────────
        run_duration = float(rt[-1] - rt[0])
        sampling_hz = float(len(rt) / run_duration) if run_duration > 0 else np.nan

        metrics = {
            "baseline_drift": round(drift_abs, 4),
            "baseline_drift_pct": round(drift_pct, 4),
            "baseline_noise_sigma": round(bl_noise, 6),
            "global_snr": round(float(global_snr), 2) if np.isfinite(global_snr) else None,
            "total_integrated_area": round(total_integrated, 2),
            "area_coverage_pct": round(area_coverage, 4),
            "n_peaks": int(len(df)),
            "mean_cqi": round(mean_cqi, 4) if np.isfinite(mean_cqi) else None,
            "overall_quality_score": round(oqs, 4),
            "run_duration_s": round(run_duration, 3),
            "sampling_rate_hz": round(sampling_hz, 3) if np.isfinite(sampling_hz) else None,
        }

        self.audit.info(
            "GlobalMetrics",
            f"Métricas globais: drift={drift_pct:.1f}%, SNR_global={global_snr:.1f}, "
            f"área_total={total_integrated:.0f}, OQS={oqs:.3f}.",
            **{k: v for k, v in metrics.items()},
        )
        return metrics

    # ==========================================================
    # 🔬  FINGERPRINTING — COMPARAÇÃO ENTRE CROMATOGRAMAS
    # ==========================================================
    @staticmethod
    def compare_chromatograms(
        rt1: np.ndarray,
        sig1: np.ndarray,
        rt2: np.ndarray,
        sig2: np.ndarray,
        *,
        interpolate: bool = True,
    ) -> dict:
        """
        Compara dois cromatogramas ponto-a-ponto (Fingerprinting).

        Esta função **não** compara listas de picos — ela compara as curvas
        completas de sinal, o que permite detectar diferenças globais (ex.:
        impurezas fora das janelas de integração, mudanças no perfil de baseline,
        eluição incompleta) que métricas por pico não capturam.

        Uso típico
        ----------
        - Amostra vs Referência: ``compare_chromatograms(rt_ref, sig_ref, rt_sample, sig_sample)``
        - Corrida A vs Corrida B (mesmo método): validar reprodutibilidade de sistema
        - Antes vs Depois de manutenção: detectar deriva instrumental

        Métricas retornadas
        -------------------
        pearson_r
            Coeficiente de correlação de Pearson entre os dois sinais
            interpolados na mesma grade de RT.
            Valores próximos de 1.0 indicam perfis idênticos.
            Sensível a deslocamentos de escala (diferenças de concentração).
        pearson_r2
            Coeficiente de determinação (R²) = pearson_r².
        rmse
            Erro médio quadrático (Root Mean Square Error).
            Mesma unidade do sinal (u.a.). Sensível a diferenças absolutas.
        nrmse
            RMSE normalizado pelo range de sig1: RMSE / (max(sig1)−min(sig1)).
            Permite comparação entre corridas de escalas diferentes.
        cosine_similarity
            Similaridade do cosseno: ``sig1 · sig2 / (‖sig1‖ × ‖sig2‖)``.
            Insensível a diferenças de escala absoluta — mede "forma" do perfil.
            1.0 = perfis identicamente proporcionais.
        spectral_contrast_angle_deg
            Ângulo entre os vetores dos dois sinais (graus).
            0° = idênticos, 90° = ortogonais (sem similaridade de forma).
            Derivado do cosine_similarity: ``arccos(cos_sim) × 180/π``.
        mae
            Mean Absolute Error — menos sensível a outliers que RMSE.
        max_abs_diff
            Diferença absoluta máxima em qualquer ponto da corrida (u.a.).
            Útil para identificar janelas de RT com maior discrepância.
        rt_max_diff
            RT onde a diferença absoluta é máxima (em segundos).
        n_points
            Número de pontos usados na comparação (após interpolação/alinhamento).
        rt_overlap_s
            Faixa de RT compartilhada entre as duas corridas (segundos).
        verdict
            Classificação qualitativa automática baseada em pearson_r:
            "IDENTICAL" (r≥0.999), "SIMILAR" (r≥0.99), "ACCEPTABLE" (r≥0.95),
            "DIFFERENT" (r<0.95).

        Parameters
        ----------
        rt1, sig1 : np.ndarray
            Tempo de retenção e sinal da corrida de referência.
        rt2, sig2 : np.ndarray
            Tempo de retenção e sinal da corrida a comparar.
        interpolate : bool
            Se True (padrão), interpola sig2 na grade de rt1
            (necessário quando as corridas têm grades de RT diferentes).
            Se False, as grades devem ter exatamente o mesmo comprimento.

        Returns
        -------
        dict
            Dicionário com todas as métricas listadas acima.

        Raises
        ------
        ValueError
            Se as corridas não tiverem sobreposição de RT suficiente (<10 pts).
        """
        from scipy.stats import pearsonr

        # ── Determina faixa de RT comum ──────────────────────────────────────
        rt_lo = max(rt1[0], rt2[0])
        rt_hi = min(rt1[-1], rt2[-1])
        rt_overlap = rt_hi - rt_lo

        if rt_overlap <= 0:
            raise ValueError(
                f"As corridas não se sobrepõem em RT " f"([{rt1[0]:.1f}–{rt1[-1]:.1f}s] vs [{rt2[0]:.1f}–{rt2[-1]:.1f}s])."
            )

        if interpolate:
            # Usa a grade de rt1 dentro da janela de sobreposição
            mask1 = (rt1 >= rt_lo) & (rt1 <= rt_hi)
            rt_common = rt1[mask1]
            a = sig1[mask1]
            b = np.interp(rt_common, rt2, sig2)
        else:
            if len(rt1) != len(rt2):
                raise ValueError(f"interpolate=False exige arrays do mesmo comprimento " f"({len(rt1)} vs {len(rt2)}).")
            rt_common = rt1
            a = sig1
            b = sig2

        if len(rt_common) < 10:
            raise ValueError(f"Sobreposição de RT insuficiente: apenas {len(rt_common)} pontos.")

        # ── Pearson ───────────────────────────────────────────────────────────
        r, _ = pearsonr(a, b)
        r2 = float(r**2)

        # ── RMSE / MAE / max_diff ─────────────────────────────────────────────
        diff = a - b
        rmse = float(np.sqrt(np.mean(diff**2)))
        mae = float(np.mean(np.abs(diff)))
        range_a = float(np.max(a) - np.min(a))
        nrmse = rmse / range_a if range_a > 0 else np.nan

        abs_diff = np.abs(diff)
        idx_max_diff = int(np.argmax(abs_diff))
        max_abs_diff = float(abs_diff[idx_max_diff])
        rt_max_diff = float(rt_common[idx_max_diff])

        # ── Cosine similarity ─────────────────────────────────────────────────
        norm_a = float(np.linalg.norm(a))
        norm_b = float(np.linalg.norm(b))
        cos_sim = float(np.dot(a, b) / (norm_a * norm_b)) if (norm_a > 0 and norm_b > 0) else np.nan
        # clip para segurança numérica antes de arccos
        cos_sim_clipped = float(np.clip(cos_sim, -1.0, 1.0))
        angle_deg = float(np.degrees(np.arccos(cos_sim_clipped))) if np.isfinite(cos_sim) else np.nan

        # ── Verdict ───────────────────────────────────────────────────────────
        if r >= 0.999:
            verdict = "IDENTICAL"
        elif r >= 0.990:
            verdict = "SIMILAR"
        elif r >= 0.950:
            verdict = "ACCEPTABLE"
        else:
            verdict = "DIFFERENT"

        return {
            "pearson_r": round(float(r), 6),
            "pearson_r2": round(r2, 6),
            "rmse": round(rmse, 4),
            "nrmse": round(float(nrmse), 6) if np.isfinite(nrmse) else None,
            "cosine_similarity": round(cos_sim, 6) if np.isfinite(cos_sim) else None,
            "spectral_contrast_angle_deg": round(angle_deg, 4) if np.isfinite(angle_deg) else None,
            "mae": round(mae, 4),
            "max_abs_diff": round(max_abs_diff, 4),
            "rt_max_diff": round(rt_max_diff, 3),
            "n_points": int(len(rt_common)),
            "rt_overlap_s": round(float(rt_overlap), 3),
            "verdict": verdict,
        }

    # ==========================================================
    # 🗒️  EXPORTAÇÃO DO AUDIT TRAIL
    # ==========================================================
    def export_audit(self, path_json: Optional[str] = None, path_csv: Optional[str] = None) -> dict:
        if path_json:
            with open(path_json, "w", encoding="utf-8") as f:
                f.write(self.audit.to_json())
            self.audit.info("Export", f"Audit trail salvo em JSON: {path_json}.")
        if path_csv:
            self.audit.to_dataframe().to_csv(path_csv, index=False, encoding="utf-8")
            self.audit.info("Export", f"Audit trail salvo em CSV: {path_csv}.")
        return {
            "run_id": self.audit.run_id,
            "method": self.method.to_dict(),
            "summary": self.audit.summary(),
            "events": self.audit.to_dict_list(),
        }

    # ==========================================================
    # 1️⃣1️⃣  PLOT
    # ==========================================================
    def plot_results(self, rt, raw, baseline, corrected, results_df):
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=rt, y=raw, mode="lines", name="Raw Signal", line=dict(width=1)))
        fig.add_trace(go.Scatter(x=rt, y=baseline, mode="lines", name="Baseline", line=dict(dash="dash")))
        fig.add_trace(go.Scatter(x=rt, y=corrected, mode="lines", name="Corrected Signal", line=dict(width=1)))
        METHOD_COLOR = {
            "EMG": "royalblue",
            "DROP_LINE": "mediumseagreen",
            "DECONVOLUTION": "mediumpurple",
            "TANGENT_SKIM_PARENT": "darkorange",
            "TANGENT_SKIM_RIDER": "tomato",
        }
        first_flags = {k: True for k in METHOD_COLOR}
        first_area = {k: True for k in METHOD_COLOR}
        first_gauss = True
        first_g_area = True
        skim_traces = getattr(self, "_skim_traces", [None] * len(results_df))
        if len(skim_traces) != len(results_df):
            skim_traces = [None] * len(results_df)
        for (idx, row), skim in zip(results_df.iterrows(), skim_traces):
            method = row.get("integration_method", "EMG")
            color = METHOD_COLOR.get(method, "gray")
            if method != "TANGENT_SKIM_RIDER" and all(k in row and pd.notna(row[k]) for k in ["A_param", "sigma", "tau"]):
                emg_curve = self.emg(rt, row["A_param"], row["rt"], row["sigma"], row["tau"])
                label = f"{method} Fit"
                fig.add_trace(
                    go.Scatter(
                        x=rt,
                        y=emg_curve,
                        mode="lines",
                        name=label if first_flags[method] else None,
                        legendgroup=label,
                        showlegend=first_flags[method],
                        line=dict(width=2, color=color),
                        opacity=0.9,
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=rt,
                        y=emg_curve,
                        mode="lines",
                        fill="tozeroy",
                        name=f"{method} Area" if first_area[method] else None,
                        legendgroup=f"{method} Area",
                        showlegend=first_area[method],
                        line=dict(width=0, color=color),
                        opacity=0.12,
                    )
                )
                first_flags[method] = False
                first_area[method] = False
            if skim is not None:
                skim_x, skim_tan = skim
                fig.add_trace(
                    go.Scatter(
                        x=skim_x,
                        y=skim_tan,
                        mode="lines",
                        name="Tangent Skim Line",
                        legendgroup="Tangent Skim Line",
                        showlegend=(idx == 0),
                        line=dict(width=2, color="gold", dash="dashdot"),
                    )
                )
            if all(k in row and pd.notna(row[k]) for k in ["gauss_A", "gauss_mu", "gauss_sigma"]):
                gauss_curve = self.gaussian(rt, row["gauss_A"], row["gauss_mu"], row["gauss_sigma"])
                fig.add_trace(
                    go.Scatter(
                        x=rt,
                        y=gauss_curve,
                        mode="lines",
                        name="Gaussian Fits" if first_gauss else None,
                        legendgroup="Gaussian Fits",
                        showlegend=first_gauss,
                        line=dict(width=2, color="tomato", dash="dot"),
                        opacity=0.9,
                    )
                )
                fig.add_trace(
                    go.Scatter(
                        x=rt,
                        y=gauss_curve,
                        mode="lines",
                        fill="tozeroy",
                        name="Gaussian Areas" if first_g_area else None,
                        legendgroup="Gaussian Areas",
                        showlegend=first_g_area,
                        line=dict(width=0, color="tomato"),
                        opacity=0.08,
                    )
                )
                first_gauss = False
                first_g_area = False
        if not results_df.empty:
            snr_text = []
            for _, row in results_df.iterrows():
                method = row.get("integration_method", "")
                snr = row.get("snr", np.nan)
                v_pct = row.get("valley_pct", np.nan)
                area = row.get("area", np.nan)
                label = f"{method}<br>SNR={snr:.1f}" if pd.notna(snr) else method
                if pd.notna(v_pct):
                    label += f"<br>%V={v_pct:.0f}%"
                if pd.notna(area):
                    label += f"<br>Area={area:.0f}"
                snr_text.append(label)
            fig.add_trace(
                go.Scatter(
                    x=results_df["marker_rt"],
                    y=results_df["marker_height"],
                    mode="markers+text",
                    name="Detected Peaks",
                    marker=dict(size=10, symbol="circle"),
                    text=snr_text,
                    textposition="top center",
                )
            )
        fig.update_layout(
            title=f"GC Analysis [{self.method.name} v{self.method.version}] — Trapezoidal Integration + Resolution Decision Tree",
            xaxis_title="Retention Time (s)",
            yaxis_title="Intensity",
            template="plotly_white",
            height=650,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        )
        fig.show()


# ══════════════════════════════════════════════════════════════════════════════
# EXEMPLO DE USO: GCAnalyzer com um único arquivo CDF (incluindo métricas estendidas e globais)
# ══════════════════════════════════════════════════════════════════════════════

# Primeiro, crie um ProcessingMethod com parâmetros personalizados, incluindo os novos para métricas estendidas.
# Cada parâmetro é explicado abaixo:
method_single = ProcessingMethod(
    name="Teste_Unico_Extendido",  # Nome do método para identificação.
    version="1.0",  # Versão do método para rastreabilidade.
    description="Método para processamento de um único arquivo CDF com métricas estendidas e globais.",  # Descrição detalhada.
    created_by="Usuario",  # Quem criou o método.
    # Parâmetros de baseline (Whittaker AsLS):
    baseline_lam=1e8,  # λ: controla a suavidade da baseline; valores maiores resultam em baseline mais suave.
    baseline_p=0.0001,  # p: controla a assimetria; valores menores forçam a baseline a ficar abaixo do sinal.
    # Estimativa de ruído:
    noise_percentile=20,  # Percentil inferior do sinal usado para estimar o ruído de referência.
    # Detecção de picos:
    snr_threshold=3.0,  # Limiar mínimo de SNR local para aceitar um pico como válido.
    min_width_seconds=1.0,  # Largura mínima de pico em segundos; filtra ruídos de alta frequência.
    min_distance_seconds=2.0,  # Distância mínima entre picos em segundos; evita detecção de picos muito próximos.
    # Decisão de sobreposição:
    rs_deconv_threshold=1.2,  # Resolução (Rs) abaixo deste valor considera os picos como sobrepostos e aplica deconvolução ou outros métodos.
    # Classificação de sobreposição:
    valley_pct_independent=85.0,  # % do vale >= este valor: picos independentes (integração separada).
    valley_pct_dropline=50.0,  # % do vale >= este valor: usa drop-line para separar picos.
    valley_pct_skim_max=25.0,  # % do vale <= este valor: candidato a tangent skim.
    height_ratio_rider=0.15,  # Razão de alturas (menor/maior) <= este valor: considera rider peak.
    # Remoção de solvente:
    solvent_rt_cutoff_s=60.0,  # Picos com RT <= este valor em segundos são candidatos a solvente e removidos.
    solvent_area_factor=5.0,  # Picos com área > este fator * mediana das áreas são removidos como solvente.
    # Alinhamento (não usado em processamento único, mas definido para completude):
    is_rt_seconds=None,  # RT esperado do padrão interno (IS) em segundos; None usa o pico de maior área automaticamente.
    is_search_window_s=10.0,  # Janela de busca ± segundos ao redor de is_rt_seconds para localizar o IS.
    rrt_bin_tolerance=0.02,  # Tolerância para agrupar picos em bins de RRT no alinhamento multi-corrida.
    # Novos parâmetros para métricas estendidas:
    dead_time_s=20.0,  # Tempo morto (t₀) da coluna em segundos; usado para k' e α (0.0 desativa).
    cqi_weight_n=1.0,  # Peso para eficiência (N_plates) no CQI.
    cqi_weight_rs=1.0,  # Peso para resolução (Rs_usp) no CQI.
    cqi_weight_tf=1.0,  # Peso para simetria (Tailing Factor) no CQI.
    cqi_weight_snr=1.0,  # Peso para SNR no CQI.
    cqi_n_ref=5000.0,  # N de referência para normalização no CQI.
    cqi_rs_ref=1.5,  # Rs mínimo aceitável para normalização no CQI.
    cqi_snr_ref=10.0,  # SNR de referência para normalização no CQI.
)

# Salve o método em JSON para reutilização (opcional, mas recomendado para rastreabilidade).
method_single.save("Metodo_Teste_Unico_Extendido.json")

# Instancie o GCAnalyzer com o método.
# Parâmetros:
# - method: O ProcessingMethod criado acima.
# - run_id: Identificador da corrida (ex.: nome do arquivo sem extensão).
# - echo_audit: Se True, imprime eventos de audit em tempo real no console.
gc_single = GCAnalyzer(method=method_single, run_id="T1", echo_audit=False)

# Caminho do arquivo CDF único.
cdf_path_single = r"C:\Users\BDanielS\Desktop\UFMG\Vault\Doutorado\Codigos\data\gc-data\Testes\cdf\T1.CDF"

# Etapa 1: Leia o CDF.
# Retorna: rt (array de tempos de retenção), intensity (array de intensidades brutas).
rt, intensity = gc_single.read_cdf(cdf_path_single)

# Etapa 2: Remova a baseline usando Whittaker AsLS (parâmetros baseline_lam e baseline_p do método).
# Retorna: corrected (sinal corrigido), baseline (array da baseline estimada).
corrected, baseline = gc_single.remove_baseline(rt, intensity)

# Etapa 3: Integre os picos (usa todos os parâmetros de detecção, sobreposição, etc. do método).
# Retorna: DataFrame com resultados dos picos (RT, área, SNR, método de integração, etc.).
results_df = gc_single.integrate(rt, corrected)

# Etapa 4: Compute métricas USP/EP como N_plates_ep, N_plates_usp, tailing_factor_usp, asymmetry_factor_ep, W_half_s, W_base_s, area_pct.
# Requer os arrays rt e corrected, e o DataFrame de resultados.
results_df = gc_single.compute_usp_metrics(rt, corrected, results_df)
print("Métricas USP/EP:")
print(results_df)  # Exibe a tabela atualizada com métricas USP/EP.

# Etapa 5: Compute métricas estendidas como Rs_usp, Rs_ep, k_prime, alpha, shape_quality_score, rsd_pct_area, CQI.
# Requer os arrays rt e intensity, o DataFrame atualizado, e opcionalmente dead_time_s (sobrescreve o do método).
# Usa dead_time_s do método se não fornecido.
results_df = gc_single.compute_extended_metrics(rt, corrected, results_df)
print("Métricas Estendidas:")
print(results_df)  # Exibe a tabela com métricas estendidas adicionadas.

# Etapa 6: Compute métricas globais da corrida como baseline_drift, global_snr, total_integrated_area, overall_quality_score (OQS).
# Requer os arrays rt, intensity (raw), corrected, baseline, e o DataFrame de resultados.
global_metrics = gc_single.compute_global_metrics(rt, intensity, corrected, baseline, results_df)
print("Métricas Globais:")
print(global_metrics)  # Exibe o dicionário com métricas globais.

# Etapa 7: Plote os resultados (sinal bruto, baseline, corrigido, fits EMG/Gauss, marcadores de picos).
gc_single.plot_results(rt, intensity, baseline, corrected, results_df)

# Etapa 8: Exporte o audit trail (log de todos os eventos e decisões).
# Salva em JSON e/ou CSV; retorna dict com o conteúdo.
audit_export = gc_single.export_audit(path_json="audit_T1_ext.json", path_csv="audit_T1_ext.csv")

# ══════════════════════════════════════════════════════════════════════════════
# EXEMPLO DE USO: GCAnalyzer com múltiplos arquivos CDF (process_batch, incluindo métricas estendidas, globais e comparação)
# ══════════════════════════════════════════════════════════════════════════════

# Crie um ProcessingMethod para o lote (pode reutilizar ou criar novo).
# Aqui, reutilizamos o mesmo, mas ajustamos parâmetros se necessário (ex.: para alinhamento preciso e métricas estendidas).
method_batch = ProcessingMethod.load("Metodo_Teste_Unico_Extendido.json")  # Carrega do JSON salvo anteriormente.
# Ajuste parâmetros específicos para lote, se quiser (ex.: defina is_rt_seconds para alinhamento preciso).
method_batch.is_rt_seconds = 82  # Exemplo: RT esperado do IS em 82s (ajuste conforme seus dados).
method_batch.rrt_bin_tolerance = 0.03  # Aumente tolerância se houver variação entre corridas.
method_batch.dead_time_s = 20.0  # Tempo morto para k' e α nas métricas estendidas.
method_batch.save("Metodo_Teste_Batch_Extendido.json")  # Salve a versão ajustada.

# Instancie o GCAnalyzer para o lote.
gc_batch = GCAnalyzer(method=method_batch, run_id="Batch_T1_T2_T3", echo_audit=False)

# Lista de arquivos CDF (use caminhos reais).
cdf_files = [
    r"C:\Users\BDanielS\Desktop\UFMG\Vault\Doutorado\Codigos\data\gc-data\Testes\cdf\T1.CDF",
    r"C:\Users\BDanielS\Desktop\UFMG\Vault\Doutorado\Codigos\data\gc-data\Testes\cdf\T2.CDF",
    r"C:\Users\BDanielS\Desktop\UFMG\Vault\Doutorado\Codigos\data\gc-data\Testes\cdf\T3.CDF",
]

# Processamento em lote:
# Parâmetros:
# - cdf_files: Lista de caminhos.
# - compute_usp: True para calcular métricas USP em cada corrida OK.
# - align: True para alinhar as corridas OK por RRT e retornar estatísticas (usa is_rt_seconds, rrt_bin_tolerance, etc.).
results_batch, df_stats = gc_batch.process_batch(cdf_files, compute_usp=True, align=True)

# Inspecione os resultados:
# results_batch é uma lista de RunResult (um por arquivo).
for res in results_batch:
    if res.ok:
        print(f"OK: {res.run_id} - {len(res.results_df)} picos")
    else:
        print(f"FAILED: {res.run_id} - {res.error_type}: {res.error_message}")

# df_stats: DataFrame de estatísticas por bin (se align=True e >=2 corridas OK).
if df_stats is not None:
    print("Estatísticas de Alinhamento:")
    print(df_stats.head())  # Exibe as primeiras linhas.
    df_stats.to_csv("stats_batch_ext.csv", index=False)  # Salve se quiser.

# Aplique métricas estendidas e globais para cada corrida OK no batch.
# Para isso, precisamos recriar os arrays rt, intensity, corrected, baseline para cada um,
# pois não são armazenados no RunResult (apenas results_df).
# Usamos uma instância temporária de GCAnalyzer para ler e remover baseline.
global_metrics_batch = {}  # Dicionário para armazenar métricas globais por run_id.
for res in results_batch:
    if res.ok:
        print(f"Calculando métricas estendidas e globais para {res.run_id}...")
        gc_temp = GCAnalyzer(method=method_batch, run_id=res.run_id, echo_audit=False)  # echo=False para não poluir console
        rt, intensity = gc_temp.read_cdf(res.cdf_path)
        corrected, baseline = gc_temp.remove_baseline(rt, intensity)

        # Métricas estendidas (requer compute_usp_metrics já chamado no batch).
        res.results_df = gc_temp.compute_extended_metrics(rt, corrected, res.results_df)
        print(f"Métricas Estendidas para {res.run_id}:")
        print(res.results_df)  # Exibe a tabela atualizada.

        # Métricas globais.
        global_metrics = gc_temp.compute_global_metrics(rt, intensity, corrected, baseline, res.results_df)
        global_metrics_batch[res.run_id] = global_metrics
        print(f"Métricas Globais para {res.run_id}:")
        print(global_metrics)

        # Plot individual.
        gc_temp.plot_results(rt, intensity, baseline, corrected, res.results_df)

# (Opcional) Compare cromatogramas entre duas corridas do batch (ex.: T1 vs T2).
# Requer recriar os arrays para as duas corridas.
if len([r for r in results_batch if r.ok]) >= 2:
    # Exemplo: Compare a primeira (ref) com a segunda.
    res1 = next(r for r in results_batch if r.ok)  # Primeira OK
    res2 = next(r for r in results_batch if r.ok and r.run_id != res1.run_id)  # Segunda OK

    gc_temp1 = GCAnalyzer(method=method_batch, run_id=res1.run_id, echo_audit=False)
    rt1, sig1 = gc_temp1.read_cdf(res1.cdf_path)  # Usa sinal bruto (intensity) para comparação.

    gc_temp2 = GCAnalyzer(method=method_batch, run_id=res2.run_id, echo_audit=False)
    rt2, sig2 = gc_temp2.read_cdf(res2.cdf_path)

    comparison = gc_temp1.compare_chromatograms(rt1, sig1, rt2, sig2, interpolate=True)
    print(f"Comparação entre {res1.run_id} e {res2.run_id}:")
    print(comparison)  # Exibe o dicionário com métricas de similaridade (pearson_r, rmse, cosine_similarity, etc.).

# Exporte o audit trail do lote (inclui eventos de todas as corridas).
audit_batch = gc_batch.export_audit(path_json="audit_batch_ext.json", path_csv="audit_batch_ext.csv")
