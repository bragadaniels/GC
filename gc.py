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
    """
    Classe base de todas as exceções do GCAnalyzer.

    Usar subclasses específicas nos blocos ``except`` permite distinguir
    falhas de arquivo (recuperável com o próximo CDF) de falhas de lógica
    (bugs reais que merecem propagar).

    Hierarquia
    ----------
    GCAnalyzerError
    ├── CDFReadError        – arquivo ausente, corrompido, variável faltante
    ├── BaselineError       – Whittaker divergiu ou produziu resultado inválido
    ├── PeakDetectionError  – nenhum pico válido encontrado após todos os filtros
    ├── IntegrationError    – falha em todas as tentativas de integração
    └── AlignmentError      – IS não encontrado, corridas incompatíveis, etc.
    """

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
# 📦  RunResult — envelope de resultado por corrida no processamento em lote
# ══════════════════════════════════════════════════════════════════════════════


@dataclass
class RunResult:
    """
    Envelope imutável que representa o resultado de uma corrida no pipeline.

    Sempre é retornado — independente de sucesso ou falha — para que o loop
    de lote nunca quebre e o chamador tenha um tipo uniforme para inspecionar.

    Atributos
    ---------
    run_id : str
        Identificador da corrida (tipicamente o nome do arquivo CDF).
    status : str
        ``"OK"`` se o pipeline completou sem erros fatais.
        ``"FAILED"`` se uma exceção não recuperável foi capturada.
    cdf_path : str
        Caminho do arquivo CDF processado.
    results_df : pd.DataFrame | None
        DataFrame de picos retornado por ``integrate()``.
        ``None`` se ``status == "FAILED"``.
    audit_events : list[dict]
        Snapshot dos eventos do audit trail desta corrida.
        Preservado mesmo em caso de falha.
    error_type : str | None
        Nome da classe de exceção capturada (ex.: ``"CDFReadError"``).
    error_message : str | None
        Mensagem da exceção.
    error_traceback : str | None
        Traceback completo como string (para log detalhado).
    """

    run_id: str
    status: str  # "OK" | "FAILED"
    cdf_path: str
    results_df: Optional[pd.DataFrame]  # None se FAILED
    audit_events: list  # sempre preenchido
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
        return f"RunResult(run_id={self.run_id!r}, status=FAILED, " f"error={self.error_type}: {self.error_message!r})"


# ══════════════════════════════════════════════════════════════════════════════
# ⚙️  PROCESSING METHOD — parâmetros versionados e rastreáveis
# ══════════════════════════════════════════════════════════════════════════════


@dataclass
class ProcessingMethod:
    """
    Encapsula todos os parâmetros do pipeline GCAnalyzer em um objeto
    versionado, auditável e serializável em JSON.

    Por quê?
    --------
    Parâmetros soltos no código (ou passados ad-hoc) tornam impossível
    garantir que um reprocessamento meses depois use *exatamente* os mesmos
    critérios. Com ProcessingMethod:

      1. Os parâmetros ficam num arquivo .json nomeado e versionado
         (ex.: ``Metodo_Benzene_V2.json``).
      2. O AuditLogger registra automaticamente qual arquivo foi usado.
      3. ``GCAnalyzer`` nunca aceita parâmetros soltos — sempre um
         ``ProcessingMethod``.

    Uso rápido
    ----------
    >>> m = ProcessingMethod(name="Benzene_V2", snr_threshold=5)
    >>> m.save("Metodo_Benzene_V2.json")
    >>> m2 = ProcessingMethod.load("Metodo_Benzene_V2.json")
    >>> analyzer = GCAnalyzer(method=m2)
    """

    # ── Metadados do método ──────────────────────────────────────────────────
    name: str = "default"
    version: str = "1.0"
    description: str = ""
    created_by: str = ""
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds"))
    source_file: str = ""  # preenchido automaticamente em load()

    # ── Baseline (Whittaker AsLS) ────────────────────────────────────────────
    baseline_lam: float = 1e8  # λ — suavização; maior = linha mais suave
    baseline_p: float = 0.0001  # p — assimetria; valores pequenos forçam baseline abaixo do sinal

    # ── Estimativa de ruído global ───────────────────────────────────────────
    noise_percentile: int = 20  # percentil inferior usado para região de referência

    # ── Detecção de picos ────────────────────────────────────────────────────
    snr_threshold: float = 3.0  # SNR local mínimo para aceitar um pico
    min_width_seconds: float = 1.0  # largura mínima (s) — evita ruído de alta frequência
    min_distance_seconds: float = 2.0  # distância mínima entre picos (s)

    # ── Decisão de sobreposição ──────────────────────────────────────────────
    rs_deconv_threshold: float = 1.2  # Rs abaixo deste valor → trata como sobreposição

    # ── Classificação de sobreposição (%Vale e razão de alturas) ────────────
    valley_pct_independent: float = 85.0  # %Vale ≥ este → picos independentes
    valley_pct_dropline: float = 50.0  # %Vale ≥ este → Drop-line
    valley_pct_skim_max: float = 25.0  # %Vale ≤ este → candidato a Tangent Skim
    height_ratio_rider: float = 0.15  # h_menor/h_maior ≤ este → Rider Peak

    # ── Remoção do pico de solvente ──────────────────────────────────────────
    solvent_rt_cutoff_s: float = 60.0  # picos com RT ≤ este valor são candidatos a solvente
    solvent_area_factor: float = 5.0  # picos com área > fator × mediana também são removidos

    # ── Alinhamento multi-corrida por RRT (Relative Retention Time) ──────────
    # RRT = RT_pico / RT_padrão_interno
    # Cancela a maior parte da deriva de fluxo e temperatura entre corridas.
    #
    # is_rt_seconds      : RT esperado do Padrão Interno em segundos.
    #                      None → usa o pico de maior área como IS automático.
    # is_search_window_s : janela ±s em torno de is_rt_seconds para localizar
    #                      o IS (evita capturar uma impureza vizinha).
    # rrt_bin_tolerance  : dois picos de corridas distintas são considerados
    #                      o mesmo composto se |RRT_a - RRT_b| ≤ este valor.
    #                      Regra prática: 0.01-0.02 para GC isotérmico estável;
    #                      até 0.05 em gradientes de temperatura ou colunas velhas.
    is_rt_seconds: Optional[float] = None
    is_search_window_s: float = 10.0
    rrt_bin_tolerance: float = 0.02

    # ────────────────────────────────────────────────────────────────────────
    # Serialização
    # ────────────────────────────────────────────────────────────────────────

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), ensure_ascii=False, indent=indent)

    def save(self, path: str | Path) -> None:
        """Salva o método em disco como JSON."""
        path = Path(path)
        path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "ProcessingMethod":
        """Carrega um método de um arquivo JSON e registra o caminho de origem."""
        path = Path(path)
        data = json.loads(path.read_text(encoding="utf-8"))
        data["source_file"] = str(path.resolve())
        return cls(**data)

    @classmethod
    def default(cls) -> "ProcessingMethod":
        """Retorna o método padrão com todos os valores nominais."""
        return cls(name="default", description="Parâmetros padrão de fábrica.")

    # ── Representação ────────────────────────────────────────────────────────

    def __str__(self) -> str:
        src = f" | fonte: {self.source_file}" if self.source_file else ""
        return f"ProcessingMethod(name={self.name!r}, version={self.version!r}" f"{src})"

    def __repr__(self) -> str:
        return self.__str__()


# ══════════════════════════════════════════════════════════════════════════════
# 🗒️  AUDIT TRAIL — estrutura de evento e logger
# ══════════════════════════════════════════════════════════════════════════════


@dataclass
class AuditEvent:
    """Representa um único evento rastreável no pipeline analítico."""

    timestamp: str
    level: str  # INFO | WARN | ERROR | DECISION
    module: str  # Baseline | PeakDetection | Integration | QC | ...
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
    """
    Audit trail estruturado para o GCAnalyzer.

    Exportação:
        logger.to_dataframe()   → pd.DataFrame
        logger.to_json()        → str JSON
        logger.to_dict_list()   → list[dict]
    """

    def __init__(self, run_id: Optional[str] = None, echo: bool = False):
        self.run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
        self.echo = echo
        self._events: list[AuditEvent] = []
        self._py_log = logging.getLogger(f"GCAnalyzer.{self.run_id}")

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat(timespec="milliseconds")

    def _add(
        self,
        level: str,
        module: str,
        message: str,
        original_value=None,
        new_value=None,
        **ctx,
    ) -> None:
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
        py_level = {
            "INFO": logging.INFO,
            "WARN": logging.WARNING,
            "ERROR": logging.ERROR,
            "DECISION": logging.INFO,
        }.get(level, logging.DEBUG)
        self._py_log.log(py_level, str(evt))

    # ── API pública ──────────────────────────────────────────────────────────

    def info(self, module: str, message: str, **ctx):
        self._add("INFO", module, message, **ctx)

    def warn(self, module: str, message: str, original_value=None, new_value=None, **ctx):
        self._add("WARN", module, message, original_value=original_value, new_value=new_value, **ctx)

    def error(self, module: str, message: str, **ctx):
        self._add("ERROR", module, message, **ctx)

    def decision(self, module: str, message: str, original_value=None, new_value=None, **ctx):
        self._add("DECISION", module, message, original_value=original_value, new_value=new_value, **ctx)

    # ── Conveniências de domínio ─────────────────────────────────────────────

    def log_method(self, method: ProcessingMethod) -> None:
        """Registra o ProcessingMethod completo como primeiro evento do pipeline."""
        src = method.source_file if method.source_file else "<in-memory>"
        self.info(
            "Method",
            f"ProcessingMethod carregado: nome='{method.name}', " f"versão='{method.version}', fonte='{src}'.",
            method_name=method.name,
            method_version=method.version,
            method_source=src,
            method_description=method.description,
            method_created_by=method.created_by,
            method_created_at=method.created_at,
            params=method.to_dict(),
        )

    def log_baseline(self, lam: float, p: float, mean_reduction: float, pts: int):
        self.info(
            "Baseline",
            f"Whittaker AsLS aplicado (λ={lam:.0e}, p={p}). "
            f"Redução média de intensidade: {mean_reduction:.1f} u.a. ({pts} pontos).",
            lam=lam,
            p=p,
            mean_reduction=mean_reduction,
            n_points=pts,
        )

    def log_noise_global(self, sigma: float, fallback: bool = False):
        tag = " [FALLBACK via diff-std]" if fallback else ""
        self.info("NoiseEstimate", f"Ruído global (MAD×1.4826){tag}: σ={sigma:.4f}", sigma=sigma)

    def log_peak_rejection(self, reason: str, rt_s: float, value: Any = None):
        self.warn(
            "PeakDetection",
            f"Pico rejeitado em RT={rt_s:.2f}s — motivo: {reason}.",
            context_rt=rt_s,
            rejection_value=value,
        )

    def log_peaks_found(self, n: int, rts: list):
        self.info(
            "PeakDetection",
            f"{n} pico(s) detectado(s) após filtragem SNR/largura/d².",
            n_peaks=n,
            retention_times_s=rts,
        )

    def log_snr_rejection(self, n_rejected: int, threshold: float):
        if n_rejected:
            self.warn(
                "PeakDetection",
                f"{n_rejected} pico(s) rejeitado(s) por SNR local < {threshold}.",
                snr_threshold=threshold,
                n_rejected=n_rejected,
            )

    def log_overlap_decision(
        self,
        rt1: float,
        rt2: float,
        Rs: float,
        valley_pct: float,
        height_ratio: float,
        method: str,
    ):
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

    def log_emg_fallback(self, rt_s: float, area_trap: float, reason: str):
        self.decision(
            "Integration",
            f"Ajuste EMG falhou em RT={rt_s:.2f}s ({reason}). " f"Revertido para integração trapezoidal (área={area_trap:.0f}).",
            original_value="EMG",
            new_value="TRAPEZOID",
            rt=rt_s,
            area_trap=area_trap,
            failure_reason=reason,
        )

    def log_integration(
        self,
        method: str,
        rt: float,
        area: float,
        snr: float,
        window: tuple,
        extra: dict = None,
    ):
        ctx = dict(rt=rt, area=area, snr=snr, window_pts=window)
        if extra:
            ctx.update(extra)
        self.info(
            "Integration",
            f"[{method}] RT={rt:.2f}s integrado: área={area:.0f}, SNR={snr:.1f}.",
            **ctx,
        )

    def log_solvent_removal(
        self,
        n_before: int,
        n_after: int,
        median_area: float,
        rt_cutoff: float,
        factor: float,
    ):
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

    def log_deconv_audit(self, area_trap: float, area_emg: float):
        pct_err = 100 * (area_emg - area_trap) / max(area_trap, 1)
        level = "WARN" if abs(pct_err) > 10 else "INFO"
        self._add(
            level,
            "Deconvolution",
            f"Área trap total={area_trap:.0f} | EMG total={area_emg:.0f} | " f"erro={pct_err:.1f}%.",
            area_trap=area_trap,
            area_emg=area_emg,
            pct_error=pct_err,
        )

    # ── Exportação ───────────────────────────────────────────────────────────

    def to_dict_list(self) -> list[dict]:
        return [e.to_dict() for e in self._events]

    def to_dataframe(self) -> pd.DataFrame:
        rows = []
        for e in self._events:
            d = e.to_dict()
            d["context"] = json.dumps(d["context"], ensure_ascii=False)
            rows.append(d)
        return pd.DataFrame(rows)

    def to_json(self, indent: int = 2) -> str:
        payload = {
            "run_id": self.run_id,
            "generated": self._now(),
            "n_events": len(self._events),
            "events": self.to_dict_list(),
        }
        return json.dumps(payload, ensure_ascii=False, indent=indent, default=str)

    def summary(self) -> dict:
        counts: dict[str, int] = {}
        for e in self._events:
            counts[e.level] = counts.get(e.level, 0) + 1
        return {"run_id": self.run_id, "total": len(self._events), **counts}

    def __len__(self) -> int:
        return len(self._events)

    def __repr__(self) -> str:
        s = self.summary()
        return (
            f"AuditLogger(run_id={s['run_id']!r}, total={s['total']}, "
            f"INFO={s.get('INFO',0)}, WARN={s.get('WARN',0)}, "
            f"ERROR={s.get('ERROR',0)}, DECISION={s.get('DECISION',0)})"
        )


# ══════════════════════════════════════════════════════════════════════════════
# 🔬  GCAnalyzer — pipeline completo com Audit Trail + ProcessingMethod
# ══════════════════════════════════════════════════════════════════════════════


from scipy.signal import find_peaks as _scipy_find_peaks, peak_widths
from scipy.optimize import curve_fit
from scipy.stats import exponnorm
from scipy.integrate import trapezoid
from pybaselines import Baseline as _PyBaseline


class GCAnalyzer:
    """
    Pipeline de análise cromatográfica GC.

    Todos os parâmetros de processamento são encapsulados em um
    ``ProcessingMethod``, garantindo rastreabilidade completa.

    Parâmetros
    ----------
    method : ProcessingMethod, optional
        Objeto com todos os parâmetros do pipeline.
        Se omitido, usa ``ProcessingMethod.default()``.
    run_id : str, optional
        Identificador único da execução (ex.: nome do CDF).
    echo_audit : bool
        Se True, imprime eventos do audit trail no stdout em tempo real.

    Exemplo
    -------
    >>> m = ProcessingMethod.load("Metodo_Benzene_V2.json")
    >>> gc = GCAnalyzer(method=m, run_id="amostra_001")
    >>> rt, raw = gc.read_cdf("amostra_001.cdf")
    >>> corrected, baseline = gc.remove_baseline(rt, raw)
    >>> results = gc.integrate(rt, corrected)
    """

    def __init__(
        self,
        method: Optional[ProcessingMethod] = None,
        run_id: Optional[str] = None,
        echo_audit: bool = False,
    ):
        self.method = method or ProcessingMethod.default()
        self.audit = AuditLogger(run_id=run_id, echo=echo_audit)

        # O primeiro evento sempre registra qual método está sendo usado.
        self.audit.log_method(self.method)

    # ── Atalhos legíveis para os parâmetros do método ────────────────────────

    @property
    def _m(self) -> ProcessingMethod:
        """Atalho interno: self._m.snr_threshold etc."""
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
            raise CDFReadError(
                f"Falha ao ler CDF '{cdf_filepath}': {exc}",
                context={"filepath": str(cdf_filepath)},
            ) from exc

        self.audit.info(
            "IO",
            f"CDF lido: {cdf_filepath} | "
            f"{len(intensity)} pontos | "
            f"RT=[{retention_time[0]:.1f}–{retention_time[-1]:.1f}]s.",
            filepath=cdf_filepath,
            n_points=len(intensity),
            rt_start=float(retention_time[0]),
            rt_end=float(retention_time[-1]),
        )
        return retention_time, intensity

    # ==========================================================
    # 2️⃣  BASELINE – WHITTAKER (AsLS)
    # ==========================================================
    def remove_baseline(self, rt, intensity):
        """Usa self.method.baseline_lam e self.method.baseline_p."""
        lam = self._m.baseline_lam
        p = self._m.baseline_p

        try:
            baseline_fitter = _PyBaseline(x_data=rt)
            baseline, _ = baseline_fitter.asls(intensity, lam=lam, p=p)
        except Exception as exc:
            raise BaselineError(
                f"Whittaker AsLS divergiu (λ={lam:.0e}, p={p}): {exc}",
                context={"lam": lam, "p": p},
            ) from exc

        if not np.all(np.isfinite(baseline)):
            raise BaselineError(
                f"Baseline contém NaN/Inf após AsLS (λ={lam:.0e}, p={p}). " "Tente aumentar λ ou diminuir p.",
                context={"lam": lam, "p": p, "n_nan": int(np.sum(~np.isfinite(baseline)))},
            )

        corrected = intensity - baseline
        corrected[corrected < 0] = 0

        mean_reduction = float(np.mean(baseline))
        self.audit.log_baseline(lam=lam, p=p, mean_reduction=mean_reduction, pts=len(intensity))
        return corrected, baseline

    # ==========================================================
    # 3️⃣  ESTIMATIVA DE RUÍDO GLOBAL (MAD)
    # ==========================================================
    def estimate_noise_level(self, intensity):
        """Usa self.method.noise_percentile."""
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
            diff = np.diff(intensity)
            sigma = np.std(diff) / np.sqrt(2) * 0.1
            fallback = True

        self.audit.log_noise_global(sigma=float(sigma), fallback=fallback)
        return sigma

    # ==========================================================
    # 3️⃣b ESTIMATIVA DE RUÍDO LOCAL (por pico)
    # ==========================================================
    def estimate_local_snr(self, rt, intensity, peak_idx, left_bound, right_bound, width_factor=3.0):
        peak_width_pts = max(right_bound - left_bound, 1)
        half_window = int(width_factor * peak_width_pts)

        l_start = max(0, left_bound - half_window)
        l_end = left_bound
        r_start = right_bound
        r_end = min(len(rt) - 1, right_bound + half_window)

        left_pts = np.arange(l_start, l_end)
        right_pts = np.arange(r_start, r_end)

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
            safe_noise = max(global_noise, 1e-9)
            self.audit.warn(
                "NoiseEstimate",
                f"Poucos pontos de referência ({len(ref_idx)}) em RT={rt[peak_idx]:.2f}s. " "Usando ruído global como fallback.",
                n_ref_pts=len(ref_idx),
                rt=float(rt[peak_idx]),
            )
            return signal / safe_noise, local_bl_val, global_noise

        x_ref = rt[ref_idx]
        y_ref = intensity[ref_idx]
        diffs = np.diff(y_ref)
        local_noise = np.std(diffs) / np.sqrt(2)

        if local_noise <= 0 or not np.isfinite(local_noise):
            local_noise = self.estimate_noise_level(intensity)
            self.audit.warn(
                "NoiseEstimate",
                f"Ruído local inválido em RT={rt[peak_idx]:.2f}s. Usando ruído global.",
                rt=float(rt[peak_idx]),
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
                "Integration",
                f"Segmento muito curto ({len(x_seg)} pts) para integração trapezoidal.",
                left=left,
                right=right,
            )
            return 0.0, x_seg, np.zeros_like(x_seg), np.zeros_like(x_seg)

        y_left = float(y_seg[0])
        y_right = float(y_seg[-1])
        baseline_virtual = np.linspace(y_left, y_right, len(x_seg))
        y_above = np.maximum(y_seg - baseline_virtual, 0.0)
        area = float(trapezoid(y_above, x_seg))

        self.audit.info(
            "Integration",
            f"Integração trapezoidal: janela=[{left}:{right}], "
            f"área={area:.2f}, y_left={y_left:.2f}, y_right={y_right:.2f}, "
            f"max_acima_bl={np.max(y_above):.2f}.",
            window=(left, right),
            area=area,
            y_left=y_left,
            y_right=y_right,
            max_above_bl=float(np.max(y_above)),
        )
        return area, x_seg, y_above, baseline_virtual

    # ==========================================================
    # 4️⃣  DETECÇÃO DE PICOS ROBUSTA
    # ==========================================================
    def find_peaks(self, rt, intensity):
        """Usa self.method.snr_threshold, min_width_seconds, min_distance_seconds."""
        snr_threshold = self._m.snr_threshold
        min_width_seconds = self._m.min_width_seconds
        min_distance_seconds = self._m.min_distance_seconds

        dt = np.mean(np.diff(rt))
        min_width_pts = max(1, int(min_width_seconds / dt))
        min_distance_pts = max(1, int(min_distance_seconds / dt))

        noise_sigma = self.estimate_noise_level(intensity)
        if noise_sigma <= 0:
            noise_sigma = np.std(intensity) * 0.01
            self.audit.warn(
                "PeakDetection",
                "noise_sigma ≤ 0 após MAD — usando 1% do desvio padrão.",
                noise_sigma=noise_sigma,
            )

        dynamic_prominence = snr_threshold * noise_sigma

        from scipy.ndimage import gaussian_filter1d

        smoothed = gaussian_filter1d(intensity, sigma=3)
        d1 = np.gradient(smoothed, rt)
        d2 = np.gradient(d1, rt)

        self.audit.info(
            "PeakDetection",
            f"Parâmetros de detecção: dt={dt:.4f}s, "
            f"min_width={min_width_pts}pts, "
            f"min_dist={min_distance_pts}pts, "
            f"prominence≥{dynamic_prominence:.2f} (SNR×σ).",
            dt=dt,
            min_width_pts=min_width_pts,
            min_distance_pts=min_distance_pts,
            prominence=dynamic_prominence,
            snr_threshold=snr_threshold,
        )

        peaks, _ = _scipy_find_peaks(
            intensity,
            prominence=dynamic_prominence,
            distance=min_distance_pts,
            width=min_width_pts,
        )

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
                self.audit.log_peak_rejection(
                    "d²<0 não satisfeito (não é mínimo de curvatura)",
                    rt_s=float(rt[p]),
                    value=float(d2[p]),
                )
                continue
            width_seconds = widths[i] * dt
            sigma_est = width_seconds / 2.355
            if sigma_est >= max_reasonable_sigma:
                self.audit.log_peak_rejection(
                    f"Largura excessiva (σ_est={sigma_est:.2f}s ≥ {max_reasonable_sigma}s)",
                    rt_s=float(rt[p]),
                    value=sigma_est,
                )
                continue
            valid.append(i)

        peaks = peaks[valid]
        left_ips = left_ips[valid]
        right_ips = right_ips[valid]
        left_base = left_base[valid]
        right_base = right_base[valid]

        snr_values = []
        for i, p in enumerate(peaks):
            snr, _, _ = self.estimate_local_snr(rt, intensity, p, left_base[i], right_base[i])
            snr_values.append(snr)
        snr_values = np.array(snr_values)

        snr_valid = snr_values >= snr_threshold
        n_rejected = int((~snr_valid).sum())
        self.audit.log_snr_rejection(n_rejected, snr_threshold)

        peaks = peaks[snr_valid]
        left_ips = left_ips[snr_valid]
        right_ips = right_ips[snr_valid]
        left_base = left_base[snr_valid]
        right_base = right_base[snr_valid]
        snr_values = snr_values[snr_valid]

        if len(peaks) == 0:
            msg = (
                f"Todos os picos foram rejeitados pelo filtro de SNR local < {snr_threshold}. "
                "Verifique snr_threshold ou a qualidade do sinal."
            )
            self.audit.warn("PeakDetection", msg, snr_threshold=snr_threshold)
            raise PeakDetectionError(msg, context={"snr_threshold": snr_threshold})

        self.audit.log_peaks_found(n=len(peaks), rts=[round(float(rt[p]), 2) for p in peaks])

        return peaks, left_ips, right_ips, snr_values, left_base, right_base

    # ==========================================================
    # 5️⃣  MODELOS — EMG + Gaussiano
    # ==========================================================
    @staticmethod
    def emg(x, A, mu, sigma, tau):
        if sigma <= 0 or tau <= 0:
            return np.zeros_like(x, dtype=float)
        K = tau / sigma
        result = A * exponnorm.pdf(x, K=K, loc=mu, scale=sigma)
        return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)

    @staticmethod
    def multi_emg(x, *params):
        y = np.zeros_like(x, dtype=float)
        for i in range(0, len(params), 4):
            A, mu, sigma, tau = params[i], params[i + 1], params[i + 2], params[i + 3]
            y += GCAnalyzer.emg(x, A, mu, sigma, tau)
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
                f"Janela muito pequena ({len(x)} pts) em RT={rt[peak_idx]:.2f}s — " "ajuste EMG abortado.",
                rt=float(rt[peak_idx]),
                window_size=len(x),
            )
            return None

        area_trap, _, y_above, bl_virtual = self.integrate_trapezoid_segment(rt, intensity, left, right)
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
                extra=dict(
                    area_emg=area_emg,
                    sigma=sigma,
                    tau=tau,
                    tailing=tau / sigma if sigma > 0 else None,
                ),
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
                # ── Metadados de localização (índices exatos no array) ────────
                # Evitam reconstrução por np.abs(rt - val).argmin() em código
                # externo, que é ambígua em regiões de coeluição.
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
            y_gauss = self.gaussian(x, Ag, mug, abs(sigmag))
            emg_params.update(
                {
                    "gauss_A": Ag,
                    "gauss_mu": mug,
                    "gauss_sigma": abs(sigmag),
                    "area_gauss": float(trapezoid(y_gauss, x)),
                }
            )
        except Exception as e:
            self.audit.warn(
                "Integration",
                f"Ajuste Gaussiano falhou em RT={rt[peak_idx]:.2f}s: {e}.",
                rt=float(rt[peak_idx]),
                reason=str(e),
            )

        return emg_params

    # ==========================================================
    # 7️⃣  MÉTRICAS DO VALE
    # ==========================================================
    def find_valley(self, intensity, peak1_idx, peak2_idx):
        region = intensity[peak1_idx:peak2_idx]
        return peak1_idx + int(np.argmin(region))

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
        """Usa os limiares definidos em self.method."""
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
            f"Drop-line iniciado: RT={rt[peak1_idx]:.2f}s e {rt[peak2_idx]:.2f}s | " f"valley_idx={valley_idx}.",
            rt1=float(rt[peak1_idx]),
            rt2=float(rt[peak2_idx]),
            valley_idx=valley_idx,
        )
        results = []

        for pk_idx, l, r in [
            (peak1_idx, left_base, valley_idx),
            (peak2_idx, valley_idx, right_base),
        ]:
            area_trap, _, _, _ = self.integrate_trapezoid_segment(rt, intensity, l, r)
            row = self.fit_emg_peak(rt, intensity, pk_idx, l, r)

            if row is None:
                x_seg = rt[l:r]
                y_seg = intensity[l:r]
                height_trap = float(np.max(y_seg)) if len(y_seg) > 0 else 0.0
                row = {
                    "A_param": np.nan,
                    "rt": rt[pk_idx],
                    "height": height_trap,
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
                    rt_s=float(rt[pk_idx]),
                    area_trap=area_trap,
                    reason="fit_emg_peak retornou None (janela drop-line)",
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
            f"Tangent Skim iniciado: pai RT={rt[parent_idx]:.2f}s | " f"rider RT={rt[rider_idx]:.2f}s.",
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
            # índice absoluto do ápice do rider no array original
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
        if front <= 0:
            return np.nan
        return float(W / (2.0 * front))

    def asymmetry_factor(self, rt, intensity, peak_idx) -> float:
        _, left_x, right_x, front, tail = self.peak_width_at_fraction(rt, intensity, peak_idx, fraction=0.10)
        if front <= 0:
            return np.nan
        return float(tail / front)

    def theoretical_plates(self, rt, intensity, peak_idx) -> float:
        W_half, *_ = self.peak_width_at_fraction(rt, intensity, peak_idx, fraction=0.50)
        t_R = float(rt[peak_idx])
        if W_half <= 0 or t_R <= 0:
            return np.nan
        return float(5.54 * (t_R / W_half) ** 2)

    def compute_usp_metrics(self, rt: np.ndarray, intensity: np.ndarray, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calcula métricas USP/EP (N, Tailing Factor, Asymmetry Factor) para
        todos os picos do DataFrame usando **os índices armazenados pelo
        pipeline**, sem nenhuma busca por RT.

        Por quê índices em vez de RT?
        ──────────────────────────────
        A função ``integrate()`` armazena ``peak_index_apex``,
        ``peak_index_start`` e ``peak_index_end`` em cada linha do DataFrame.
        Usar esses valores elimina dois problemas:

        1. **Ineficiência** — evita refazer ``np.abs(rt - val).argmin()``
           para cada pico, operação O(N) desnecessária.
        2. **Ambiguidade** — em regiões de coeluição severa, dois picos podem
           estar tão próximos que a busca pelo RT mais próximo retorna o índice
           errado. O índice guardado é inequívoco.

        Parameters
        ----------
        rt : np.ndarray
            Array de tempo de retenção (mesma instância usada em ``integrate``).
        intensity : np.ndarray
            Sinal corrigido de baseline (mesma instância usada em ``integrate``).
        df : pd.DataFrame
            DataFrame retornado por ``integrate()``.
            Deve conter ``peak_index_apex``, ``peak_index_start``,
            ``peak_index_end``.

        Returns
        -------
        pd.DataFrame
            Cópia do DataFrame de entrada com as colunas adicionadas:
            ``N_plates``, ``tailing_factor_usp``, ``asymmetry_factor_ep``.

        Raises
        ------
        KeyError
            Se o DataFrame não contiver as colunas de índice necessárias.
        """
        required = {"peak_index_apex", "peak_index_start", "peak_index_end"}
        missing = required - set(df.columns)
        if missing:
            raise KeyError(
                f"DataFrame não contém as colunas de índice: {missing}. "
                "Certifique-se de usar o DataFrame retornado por integrate()."
            )

        out = df.copy()
        n_plates_col, tf_col, af_col = [], [], []

        for _, row in out.iterrows():
            apex = int(row["peak_index_apex"])

            if apex < 0 or apex >= len(rt):
                self.audit.warn(
                    "USPMetrics",
                    f"peak_index_apex={apex} fora do intervalo do array " f"(len={len(rt)}). Métricas definidas como NaN.",
                    apex=apex,
                )
                n_plates_col.append(np.nan)
                tf_col.append(np.nan)
                af_col.append(np.nan)
                continue

            n = self.theoretical_plates(rt, intensity, apex)
            tf = self.tailing_factor(rt, intensity, apex)
            af = self.asymmetry_factor(rt, intensity, apex)

            n_plates_col.append(n)
            tf_col.append(tf)
            af_col.append(af)

            self.audit.info(
                "USPMetrics",
                f"RT={rt[apex]:.2f}s (idx={apex}): " f"N={n:.0f}, TF={tf:.3f}, AF={af:.3f}.",
                rt=float(rt[apex]),
                apex_idx=apex,
                N=n,
                tailing_factor=tf,
                asymmetry_factor=af,
            )

        out["N_plates"] = n_plates_col
        out["tailing_factor_usp"] = tf_col
        out["asymmetry_factor_ep"] = af_col
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

        area_total_trap, _, y_above_total, _ = self.integrate_trapezoid_segment(rt, intensity, left, right)

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
            self.audit.error("Deconvolution", f"Ajuste Multi-EMG falhou: {e}. Deconvolução abortada.", reason=str(e))
            return None

        popt_g = None
        try:
            popt_g, _ = curve_fit(
                self.multi_gaussian, x, y, p0=initial_guess_g, bounds=(lower_bounds_g, upper_bounds_g), maxfev=10000
            )
        except Exception as e:
            self.audit.warn("Deconvolution", f"Ajuste Multi-Gaussiano falhou: {e}. Continuando apenas com EMG.", reason=str(e))

        emg_areas = []
        for i in range(num_peaks):
            A, mu, sigma, tau = (popt[i * 4], popt[i * 4 + 1], abs(popt[i * 4 + 2]), abs(popt[i * 4 + 3]))
            y_comp = self.emg(x, A, mu, sigma, tau)
            emg_areas.append(float(trapezoid(y_comp, x)))

        total_emg_area = sum(emg_areas)
        self.audit.log_deconv_audit(area_total_trap, total_emg_area)

        results = []
        for i in range(num_peaks):
            A, mu, sigma, tau = (popt[i * 4], popt[i * 4 + 1], abs(popt[i * 4 + 2]), abs(popt[i * 4 + 3]))
            y_comp = self.emg(x, A, mu, sigma, tau)

            frac_emg = emg_areas[i] / total_emg_area if total_emg_area > 0 else 1.0 / num_peaks
            area_comp_trap = frac_emg * area_total_trap

            self.audit.info(
                "Deconvolution",
                f"Componente {i+1}: RT={mu:.2f}s, fração EMG={frac_emg:.3f}, " f"área proporcional={area_comp_trap:.0f}.",
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
                # O índice do ápice é o pico detectado pelo scipy.
                # start/end cobrem a janela global compartilhada pelos dois picos.
                "peak_index_apex": int(peak_indices[i]),
                "peak_index_start": int(left),
                "peak_index_end": int(right),
            }
            if popt_g is not None:
                row.update(
                    {
                        "gauss_A": popt_g[i * 3],
                        "gauss_mu": popt_g[i * 3 + 1],
                        "gauss_sigma": abs(popt_g[i * 3 + 2]),
                    }
                )
            results.append(row)

        return results

    # ==========================================================
    # 9️⃣  REMOÇÃO DE PICO DO SOLVENTE
    # ==========================================================
    def remove_solvent_peak(self, df):
        """Usa self.method.solvent_rt_cutoff_s e solvent_area_factor."""
        if df.empty:
            return df
        rt_min_exclude = self._m.solvent_rt_cutoff_s
        area_factor = self._m.solvent_area_factor
        median_area = np.median(df["area"])
        filtered = df[(df["rt"] > rt_min_exclude) & (df["area"] < area_factor * median_area)]
        self.audit.log_solvent_removal(
            n_before=len(df),
            n_after=len(filtered),
            median_area=float(median_area),
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

                    else:  # DECONVOLUTION
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
                            f"{len(overlap)} pico(s) adicionado(s) via " f"{method} (RT≈{t1:.1f}–{t2:.1f}s).",
                            method=method,
                            n_peaks=len(overlap),
                            rt1=float(t1),
                            rt2=float(t2),
                        )
                        i += 2
                        continue

            # ── Pico isolado ──────────────────────────────────────────────
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
                    f"Pico isolado integrado: RT={peak_result['rt']:.2f}s, " f"área={peak_result['area']:.0f}, SNR={snr:.1f}.",
                    rt=peak_result["rt"],
                    area=peak_result["area"],
                    snr=snr,
                )
            i += 1

        # Extrai skim traces internos
        self._skim_traces = []
        clean_results = []
        for row in results:
            skim_x = row.pop("_skim_x", None)
            skim_tan = row.pop("_skim_tangent", None)
            self._skim_traces.append((skim_x, skim_tan) if skim_x is not None else None)
            clean_results.append(row)

        df = pd.DataFrame(clean_results)
        df = self.remove_solvent_peak(df)

        self.audit.info(
            "Pipeline",
            f"Integração concluída: {len(df)} pico(s) no relatório final.",
            n_peaks_final=len(df),
        )
        return df

    # ==========================================================
    # 🔗  ALINHAMENTO MULTI-CORRIDA — RRT + BINNING
    # ==========================================================

    def find_internal_standard(self, df: pd.DataFrame) -> pd.Series:
        """
        Localiza o pico do Padrão Interno (IS) em um DataFrame de corrida única.

        Estratégia de busca (em ordem de prioridade):
        ─────────────────────────────────────────────
        1. Se ``method.is_rt_seconds`` estiver definido:
           Busca o pico de **maior área** dentro da janela
           [is_rt_seconds ± is_search_window_s].
           → Robusto a pequenas flutuações de RT entre corridas.
        2. Fallback automático (``is_rt_seconds is None``):
           Usa o pico de maior área em toda a corrida.
           → Adequado para desenvolver métodos ou triagem inicial.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame retornado por ``integrate()``.
            Deve conter as colunas ``rt`` e ``area``.

        Returns
        -------
        pd.Series
            Linha do DataFrame correspondente ao IS encontrado.

        Raises
        ------
        ValueError
            Se nenhum pico for encontrado na janela especificada.
        """
        if df.empty:
            raise ValueError("DataFrame vazio — não é possível localizar o IS.")

        m = self._m

        if m.is_rt_seconds is not None:
            lo = m.is_rt_seconds - m.is_search_window_s
            hi = m.is_rt_seconds + m.is_search_window_s
            window = df[(df["rt"] >= lo) & (df["rt"] <= hi)]

            if window.empty:
                raise ValueError(
                    f"Nenhum pico encontrado na janela IS "
                    f"[{lo:.1f}–{hi:.1f}]s. "
                    f"Verifique is_rt_seconds ({m.is_rt_seconds}s) e "
                    f"is_search_window_s ({m.is_search_window_s}s)."
                )
            is_row = window.loc[window["area"].idxmax()]
            strategy = f"maior área em [{lo:.1f}–{hi:.1f}]s"
        else:
            is_row = df.loc[df["area"].idxmax()]
            strategy = "maior área global (fallback automático)"

        self.audit.info(
            "RRTAlignment",
            f"IS localizado: RT={is_row['rt']:.3f}s, área={is_row['area']:.0f} " f"({strategy}).",
            is_rt=float(is_row["rt"]),
            is_area=float(is_row["area"]),
            strategy=strategy,
        )
        return is_row

    # ----------------------------------------------------------

    def compute_rrt(self, df: pd.DataFrame, is_rt: float) -> pd.DataFrame:
        """
        Adiciona a coluna ``rrt`` (Relative Retention Time) ao DataFrame.

        RRT = RT_pico / RT_is

        A normalização pelo IS cancela as flutuações sistemáticas de fluxo
        e temperatura que deslocam **todos** os picos de uma corrida na
        mesma direção. O resultado é uma escala adimensional estável o
        suficiente para comparar corridas separadas por dias ou semanas.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame de uma única corrida (saída de ``integrate()``).
        is_rt : float
            Tempo de retenção do IS **nesta corrida** (em segundos).
            Obtido tipicamente via ``find_internal_standard(df)["rt"]``.

        Returns
        -------
        pd.DataFrame
            Cópia com a coluna ``rrt`` adicionada.

        Raises
        ------
        ValueError
            Se ``is_rt`` for zero ou negativo.
        """
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

    # ----------------------------------------------------------

    def align_runs(
        self,
        runs: list[tuple[str, pd.DataFrame]],
        auto_is: bool = True,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Alinha múltiplas corridas cromatográficas por Binning de RRT.

        O algoritmo em 4 passos
        ───────────────────────
        1. **Localizar o IS** em cada corrida (via ``find_internal_standard``).
        2. **Calcular o RRT** de todos os picos de cada corrida (via ``compute_rrt``).
        3. **Binning guloso por RRT**: percorre os picos em ordem crescente de
           RRT e agrupa num mesmo *bin* todos os picos cujo RRT está a menos de
           ``method.rrt_bin_tolerance`` do centróide atual do bin.
           O centróide é atualizado incrementalmente a cada pico adicionado.
        4. **Construir as tabelas de saída**: tabela longa (um registro por
           corrida × bin) e tabela de estatísticas (média, CV%, n, status).

        Parameters
        ----------
        runs : list of (run_id, df)
            Cada tupla é o identificador textual da corrida e o DataFrame
            retornado por ``integrate()``.
            Mínimo: 1 corrida (tabela longa útil; estatísticas triviais).
        auto_is : bool
            Se True (padrão), chama ``find_internal_standard`` em cada corrida
            para determinar o IS automaticamente conforme ``method``.
            Se False, o chamador deve garantir que os DataFrames já contenham
            a coluna ``rrt``.

        Returns
        -------
        df_long : pd.DataFrame
            Tabela longa: uma linha por (bin_id, run_id, pico).
            Colunas principais: ``bin_id``, ``bin_rrt_centroid``,
            ``run_id``, ``rt``, ``rrt``, ``area``, ``integration_method``.

        df_stats : pd.DataFrame
            Estatísticas por bin: uma linha por bin.
            Colunas: ``bin_id``, ``bin_rrt_centroid``, ``n_runs_present``,
            ``rrt_mean``, ``rrt_cv_pct``,
            ``area_mean``, ``area_cv_pct``, ``area_status``,
            ``rt_mean_s``, ``rt_std_s``.

        Raises
        ------
        ValueError
            Se ``runs`` estiver vazio.

        Notes
        -----
        *Por que binning guloso e não clustering k-means?*
        O binning guloso é determinístico, sem hiperparâmetro de número de
        clusters, e preserva a ordem cromatográfica natural dos picos.
        Para misturas com muitos componentes próximos, diminua
        ``rrt_bin_tolerance``; para colunas muito instáveis, aumente-o.
        """
        if not runs:
            raise ValueError("Lista de corridas vazia.")

        m = self._m
        tol = m.rrt_bin_tolerance

        self.audit.info(
            "RRTAlignment",
            f"align_runs iniciado: {len(runs)} corridas, " f"tolerância RRT={tol}.",
            n_runs=len(runs),
            rrt_bin_tolerance=tol,
            auto_is=auto_is,
        )

        # ── Passo 1 e 2: IS + RRT por corrida ───────────────────────────────
        enriched: list[pd.DataFrame] = []
        for run_id, df in runs:
            if auto_is:
                try:
                    is_row = self.find_internal_standard(df)
                    is_rt = float(is_row["rt"])
                except ValueError as exc:
                    self.audit.error(
                        "RRTAlignment",
                        f"Corrida '{run_id}': IS não encontrado — {exc}. " "Corrida ignorada no alinhamento.",
                        run_id=run_id,
                        error=str(exc),
                    )
                    continue
                rrt_df = self.compute_rrt(df, is_rt)
            else:
                if "rrt" not in df.columns:
                    self.audit.error(
                        "RRTAlignment",
                        f"Corrida '{run_id}' não contém coluna 'rrt' e " "auto_is=False. Corrida ignorada.",
                        run_id=run_id,
                    )
                    continue
                rrt_df = df.copy()

            rrt_df["run_id"] = run_id
            enriched.append(rrt_df)

        if not enriched:
            raise ValueError("Nenhuma corrida pôde ser processada. " "Verifique os parâmetros do IS e os DataFrames de entrada.")

        # ── Passo 3: Binning guloso por RRT ─────────────────────────────────
        # Reúne todos os picos de todas as corridas num único DataFrame
        # ordenado por RRT crescente.
        all_peaks = pd.concat(enriched, ignore_index=True).sort_values("rrt").reset_index(drop=True)

        bin_ids: list[int] = []
        centroids: list[float] = []

        current_bin = -1
        current_centroid = -np.inf
        current_count = 0

        for rrt_val in all_peaks["rrt"]:
            if abs(rrt_val - current_centroid) <= tol:
                # Pico pertence ao bin atual → atualiza centróide
                current_centroid = (current_centroid * current_count + rrt_val) / (current_count + 1)
                current_count += 1
            else:
                # Pico está além da tolerância → abre novo bin
                current_bin += 1
                current_centroid = float(rrt_val)
                current_count = 1

            bin_ids.append(current_bin)
            centroids.append(current_centroid)

        all_peaks["bin_id"] = bin_ids
        all_peaks["bin_rrt_centroid"] = [
            # centróide final do bin (após todos os picos serem alocados)
            all_peaks.loc[all_peaks["bin_id"] == b, "rrt"].mean()
            for b in bin_ids
        ]

        n_bins = all_peaks["bin_id"].nunique()
        self.audit.info(
            "RRTAlignment",
            f"Binning concluído: {len(all_peaks)} picos agrupados em " f"{n_bins} bins (tol={tol}).",
            n_peaks_total=len(all_peaks),
            n_bins=n_bins,
        )

        # ── Passo 4a: Tabela longa ───────────────────────────────────────────
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
        cols_long_present = [c for c in cols_long if c in all_peaks.columns]
        df_long = all_peaks[cols_long_present].copy()

        # ── Passo 4b: Estatísticas por bin ──────────────────────────────────
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

            rt_mean_s = float(group["rt"].mean())
            rt_std_s = float(group["rt"].std(ddof=1)) if len(group) > 1 else 0.0

            self.audit.info(
                "RRTAlignment",
                f"Bin {bin_id}: centróide RRT={centroid:.4f}, "
                f"presente em {n_present}/{n_runs_total} corridas, "
                f"CV(área)={area_cv:.1f}%.",
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
                    "rt_mean_s": rt_mean_s,
                    "rt_std_s": rt_std_s,
                }
            )

        df_stats = pd.DataFrame(stat_rows)

        self.audit.info(
            "RRTAlignment",
            f"align_runs concluído. "
            f"Bins totais: {n_bins} | "
            f"Bins presentes em todas as corridas: "
            f"{int((df_stats['n_runs_present'] == n_runs_total).sum())}.",
            n_bins=n_bins,
            n_bins_complete=int((df_stats["n_runs_present"] == n_runs_total).sum()),
        )

        return df_long, df_stats

    # ==========================================================
    # 🔄  PROCESSAMENTO EM LOTE — resiliente a falhas por arquivo
    # ==========================================================

    def process_batch(
        self,
        cdf_files: list[str | Path],
        *,
        compute_usp: bool = True,
        align: bool = False,
    ) -> tuple[list[RunResult], pd.DataFrame | None]:
        """
        Processa uma lista de arquivos CDF de forma resiliente.

        Se um arquivo falhar **em qualquer estágio**, o erro é:
        1. Registrado no Audit Trail com nível ERROR e traceback completo.
        2. Encapsulado em um ``RunResult(status="FAILED")``.
        3. A memória do objeto parcial é liberada via ``gc.collect()``.
        4. O loop avança imediatamente para o próximo arquivo.

        Os arquivos subsequentes **nunca são afetados** pela falha de um arquivo
        anterior — mesmo que o erro seja ``MemoryError`` ou corrupção de CDF.

        Hierarquia de captura por estágio
        ───────────────────────────────────
        ``CDFReadError``      → arquivo ausente/corrompido
        ``BaselineError``     → Whittaker divergiu ou produziu NaN
        ``PeakDetectionError``→ nenhum pico sobreviveu aos filtros
        ``IntegrationError``  → falha total na integração
        ``MemoryError``       → sem RAM; libera e continua
        ``GCAnalyzerError``   → qualquer erro interno tipado
        ``Exception``         → guarda-chuva final (bugs inesperados)

        Parameters
        ----------
        cdf_files : list of str or Path
            Caminhos dos arquivos CDF a processar (ordem preservada).
        compute_usp : bool
            Se True (padrão), chama ``compute_usp_metrics`` em cada corrida
            bem-sucedida.
        align : bool
            Se True, chama ``align_runs`` ao final com as corridas OK e
            retorna ``(results, df_stats)`` como segundo elemento da tupla.
            Se False, o segundo elemento é ``None``.

        Returns
        -------
        results : list[RunResult]
            Um ``RunResult`` por arquivo CDF, na mesma ordem de ``cdf_files``.
            Inspecione ``result.ok`` para distinguir sucesso de falha.
        df_stats : pd.DataFrame | None
            Tabela de estatísticas de alinhamento (apenas se ``align=True``
            e pelo menos duas corridas OK existirem). Caso contrário, ``None``.

        Examples
        --------
        >>> results, df_stats = gc.process_batch(glob("*.cdf"), align=True)
        >>> ok  = [r for r in results if r.ok]
        >>> bad = [r for r in results if not r.ok]
        >>> print(f"{len(ok)} OK | {len(bad)} falhas")
        """
        n_total = len(cdf_files)
        self.audit.info(
            "Batch",
            f"Lote iniciado: {n_total} arquivo(s) CDF.",
            n_files=n_total,
            compute_usp=compute_usp,
            align=align,
        )

        results: list[RunResult] = []

        for file_idx, cdf_path in enumerate(cdf_files, start=1):
            cdf_path = Path(cdf_path)
            run_id = cdf_path.stem
            self.audit.info(
                "Batch",
                f"[{file_idx}/{n_total}] Iniciando: {cdf_path.name}",
                run_id=run_id,
                file_idx=file_idx,
            )

            # Variáveis locais — limpas em caso de falha
            rt = intensity = corrected = baseline = df = None

            try:
                # ── Leitura ──────────────────────────────────────────────────
                rt, intensity = self.read_cdf(str(cdf_path))

                # ── Baseline ─────────────────────────────────────────────────
                corrected, baseline = self.remove_baseline(rt, intensity)

                # ── Integração ───────────────────────────────────────────────
                df = self.integrate(rt, corrected)

                if df.empty:
                    raise IntegrationError(
                        "integrate() retornou DataFrame vazio — nenhum pico " "sobreviveu à remoção do solvente.",
                        context={"run_id": run_id},
                    )

                # ── Métricas USP (opcional) ───────────────────────────────────
                if compute_usp:
                    df = self.compute_usp_metrics(rt, corrected, df)

                self.audit.info(
                    "Batch",
                    f"[{file_idx}/{n_total}] OK: {len(df)} pico(s) — {cdf_path.name}",
                    run_id=run_id,
                    n_peaks=len(df),
                )

                results.append(
                    RunResult(
                        run_id=run_id,
                        status="OK",
                        cdf_path=str(cdf_path),
                        results_df=df,
                        audit_events=self.audit.to_dict_list(),
                    )
                )

            # ── Captura hierárquica — do mais específico ao mais geral ───────
            except CDFReadError as exc:
                self._handle_batch_error(
                    results,
                    run_id,
                    cdf_path,
                    exc,
                    file_idx,
                    n_total,
                    "Arquivo CDF ilegível — pulando para o próximo.",
                )

            except BaselineError as exc:
                self._handle_batch_error(
                    results,
                    run_id,
                    cdf_path,
                    exc,
                    file_idx,
                    n_total,
                    "Falha na subtração de baseline.",
                )

            except PeakDetectionError as exc:
                self._handle_batch_error(
                    results,
                    run_id,
                    cdf_path,
                    exc,
                    file_idx,
                    n_total,
                    "Nenhum pico detectado.",
                )

            except IntegrationError as exc:
                self._handle_batch_error(
                    results,
                    run_id,
                    cdf_path,
                    exc,
                    file_idx,
                    n_total,
                    "Falha na integração.",
                )

            except GCAnalyzerError as exc:
                self._handle_batch_error(
                    results,
                    run_id,
                    cdf_path,
                    exc,
                    file_idx,
                    n_total,
                    "Erro interno do GCAnalyzer.",
                )

            except MemoryError as exc:
                # MemoryError não herda de GCAnalyzerError — tratamento especial
                self._handle_batch_error(
                    results,
                    run_id,
                    cdf_path,
                    exc,
                    file_idx,
                    n_total,
                    "Memória insuficiente para processar este arquivo.",
                )

            except Exception as exc:
                # Guarda-chuva final — captura bugs inesperados sem derrubar o lote
                self._handle_batch_error(
                    results,
                    run_id,
                    cdf_path,
                    exc,
                    file_idx,
                    n_total,
                    "Erro inesperado (verifique o traceback no audit trail).",
                )

            finally:
                # Libera memória dos arrays deste arquivo independente do resultado
                del rt, intensity, corrected, baseline, df
                _gc.collect()

        n_ok = sum(1 for r in results if r.ok)
        n_fail = len(results) - n_ok
        self.audit.info(
            "Batch",
            f"Lote concluído: {n_ok} OK | {n_fail} falha(s) | " f"{n_total} total.",
            n_ok=n_ok,
            n_fail=n_fail,
            n_total=n_total,
        )

        # ── Alinhamento opcional ──────────────────────────────────────────────
        df_stats: pd.DataFrame | None = None
        if align:
            ok_runs = [(r.run_id, r.results_df) for r in results if r.ok]
            if len(ok_runs) >= 2:
                try:
                    _, df_stats = self.align_runs(ok_runs)
                    self.audit.info(
                        "Batch",
                        f"Alinhamento concluído sobre {len(ok_runs)} corridas OK.",
                        n_aligned=len(ok_runs),
                    )
                except AlignmentError as exc:
                    self.audit.error(
                        "Batch",
                        f"Alinhamento falhou após lote: {exc}.",
                        error=str(exc),
                    )
            else:
                self.audit.warn(
                    "Batch",
                    f"Alinhamento solicitado mas apenas {len(ok_runs)} corrida(s) OK " f"(mínimo: 2). Ignorado.",
                    n_ok_runs=len(ok_runs),
                )

        return results, df_stats

    # ----------------------------------------------------------

    def _handle_batch_error(
        self,
        results: list,
        run_id: str,
        cdf_path: Path,
        exc: BaseException,
        file_idx: int,
        n_total: int,
        summary: str,
    ) -> None:
        """
        Centraliza o tratamento de erro de uma corrida no lote:
        loga no audit trail (ERROR), anexa o RunResult FAILED, e retorna.
        Nunca re-levanta — garante que o loop continue.
        """
        tb_str = traceback.format_exc()
        exc_type = type(exc).__name__

        self.audit.error(
            "Batch",
            f"[{file_idx}/{n_total}] FAILED — {run_id}: {summary} " f"({exc_type}: {exc})",
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
            "method": self.method.to_dict(),  # ← método snapshottado junto
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
            title=(
                f"GC Analysis [{self.method.name} v{self.method.version}] — " "Trapezoidal Integration + Resolution Decision Tree"
            ),
            xaxis_title="Retention Time (s)",
            yaxis_title="Intensity",
            template="plotly_white",
            height=650,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        )
        fig.show()


# ══════════════════════════════════════════════════════════════════════════════
# EXEMPLO DE USO: GCAnalyzer com um único arquivo CDF
# ══════════════════════════════════════════════════════════════════════════════

# Primeiro, crie um ProcessingMethod com parâmetros personalizados.
# Cada parâmetro é explicado abaixo:
method_single = ProcessingMethod(
    name="Teste_Unico",  # Nome do método para identificação.
    version="1.0",  # Versão do método para rastreabilidade.
    description="Método para processamento de um único arquivo CDF de teste.",  # Descrição detalhada.
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
)

# Salve o método em JSON para reutilização (opcional, mas recomendado para rastreabilidade).
method_single.save("Metodo_Teste_Unico.json")

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

# (Opcional) Compute métricas USP/EP como N_plates, tailing_factor_usp, asymmetry_factor_ep.
# Requer os arrays rt e corrected, e o DataFrame de resultados.
results_df = gc_single.compute_usp_metrics(rt, corrected, results_df)
print(results_df)  # Exibe a tabela atualizada com métricas USP/EP.

# Etapa 4: Plote os resultados (sinal bruto, baseline, corrigido, fits EMG/Gauss, marcadores de picos).
gc_single.plot_results(rt, intensity, baseline, corrected, results_df)

# Etapa 5: Exporte o audit trail (log de todos os eventos e decisões).
# Salva em JSON e/ou CSV; retorna dict com o conteúdo.
audit_export = gc_single.export_audit(path_json="audit_T1.json", path_csv="audit_T1.csv")

# ══════════════════════════════════════════════════════════════════════════════
# EXEMPLO DE USO: GCAnalyzer com múltiplos arquivos CDF (process_batch)
# ══════════════════════════════════════════════════════════════════════════════

# Crie um ProcessingMethod para o lote (pode reutilizar ou criar novo).
# Aqui, reutilizamos o mesmo, mas ajustamos parâmetros se necessário (ex.: para alinhamento).
method_batch = ProcessingMethod.load("Metodo_Teste_Unico.json")  # Carrega do JSON salvo anteriormente.
# Ajuste parâmetros específicos para lote, se quiser (ex.: defina is_rt_seconds para alinhamento preciso).
method_batch.is_rt_seconds = None  # Exemplo: RT esperado do IS em 120s (ajuste conforme seus dados).
method_batch.rrt_bin_tolerance = 0.03  # Aumente tolerância se houver variação entre corridas.
method_batch.save("Metodo_Teste_Batch.json")  # Salve a versão ajustada.

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
    df_stats.to_csv("stats_batch.csv", index=False)  # Salve se quiser.

# Plot individual para cada corrida OK no batch.
# Para plotar, precisamos recriar os arrays rt, intensity, corrected, baseline para cada um,
# pois não são armazenados no RunResult (apenas results_df).
# Usamos uma instância temporária de GCAnalyzer para ler e remover baseline.
for res in results_batch:
    if res.ok:
        print(f"Plotando {res.run_id}...")
        gc_plot = GCAnalyzer(method=method_batch, run_id=res.run_id, echo_audit=False)  # echo=False para não poluir console
        rt, intensity = gc_plot.read_cdf(res.cdf_path)
        corrected, baseline = gc_plot.remove_baseline(rt, intensity)
        gc_plot.plot_results(rt, intensity, baseline, corrected, res.results_df)

# Exporte o audit trail do lote (inclui eventos de todas as corridas).
audit_batch = gc_batch.export_audit(path_json="audit_batch.json", path_csv="audit_batch.csv")
