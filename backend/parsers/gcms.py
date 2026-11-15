import io
import numpy as np
from scipy.io import netcdf_file

# Variáveis do padrão ANDI-MS que realmente usamos
_NEEDED = [
    "scan_acquisition_time",
    "total_intensity",
    "scan_index",
    "point_count",
    "mass_values",
    "intensity_values",
]


def _attr(var, name, default):
    """Lê scale_factor / add_offset de uma variável NetCDF com segurança."""
    val = getattr(var, name, None)
    if val is None and hasattr(var, "_attributes"):
        val = var._attributes.get(name)
    if val is None:
        return default
    return float(np.asarray(val).reshape(-1)[0])


def _extract_scipy(file_bytes: bytes) -> dict:
    nc = netcdf_file(io.BytesIO(file_bytes), "r", mmap=False)
    try:
        out = {}
        for name in _NEEDED:
            if name in nc.variables:
                var = nc.variables[name]
                out[name] = (
                    np.array(var[:]),
                    _attr(var, "scale_factor", 1.0),
                    _attr(var, "add_offset", 0.0),
                )
        return out
    finally:
        nc.close()


def _extract_netcdf4(file_bytes: bytes, filename: str) -> dict:
    import netCDF4

    nc = netCDF4.Dataset(filename, memory=file_bytes)
    try:
        # Desliga o auto-scale: a escala é aplicada manualmente em um único lugar
        nc.set_auto_maskandscale(False)
        out = {}
        for name in _NEEDED:
            if name in nc.variables:
                var = nc.variables[name]
                out[name] = (
                    np.array(var[:]),
                    _attr(var, "scale_factor", 1.0),
                    _attr(var, "add_offset", 0.0),
                )
        return out
    finally:
        nc.close()


def read_gcms_cdf(file_bytes: bytes, filename: str) -> dict:
    """
    Lê um .CDF de GC-MS (ANDI-MS).

    Estrutura do formato: os espectros de todos os scans ficam concatenados em
    `mass_values` / `intensity_values`; `scan_index` é a posição inicial de cada scan
    e `point_count` o número de pontos (m/z) daquele scan.

    Retorna o cromatograma de íons totais (TIC) para o gráfico e, em "_arrays",
    os dados completos para servir os espectros sob demanda.
    """
    try:
        raw = _extract_scipy(file_bytes)
    except Exception:
        raw = _extract_netcdf4(file_bytes, filename)

    def get(name):
        if name not in raw:
            return None
        arr, scale, offset = raw[name]
        arr = np.asarray(arr, dtype=np.float64).reshape(-1)
        if scale != 1.0 or offset != 0.0:
            arr = arr * scale + offset
        return arr

    for required in ("scan_acquisition_time", "scan_index", "mass_values", "intensity_values"):
        if required not in raw:
            raise ValueError(f"Variável '{required}' ausente — o arquivo não parece ser um GC-MS (ANDI-MS).")

    t_sec = get("scan_acquisition_time")
    mass = get("mass_values")
    inten = get("intensity_values")
    start = get("scan_index").astype(np.int64)

    n_scans = min(len(t_sec), len(start))
    t_sec, start = t_sec[:n_scans], start[:n_scans]
    if n_scans == 0:
        raise ValueError("O arquivo não contém scans.")

    n_pts = len(mass)
    pc = get("point_count")
    if pc is not None and len(pc) >= n_scans:
        count = pc[:n_scans].astype(np.int64)
    else:
        count = np.append(start[1:], n_pts) - start
    count = np.maximum(count, 0)

    # Garante que nenhum scan aponte para fora dos arrays
    start = np.clip(start, 0, n_pts)
    count = np.minimum(count, n_pts - start)

    # TIC calculado a partir dos espectros (fallback se o arquivo não trouxer total_intensity)
    cs = np.concatenate(([0.0], np.cumsum(inten)))
    tic_calc = cs[start + count] - cs[start]
    tic = get("total_intensity")
    if tic is None or len(tic) < n_scans:
        tic = tic_calc
    else:
        tic = tic[:n_scans]

    time_min = t_sec / 60.0

    return {
        "filename": filename,
        "n_scans": int(n_scans),
        "time": np.round(time_min, 5).tolist(),
        "tic": np.round(tic, 2).tolist(),
        "mz_min": float(np.min(mass)) if n_pts else 0.0,
        "mz_max": float(np.max(mass)) if n_pts else 0.0,
        "_arrays": {
            "time0": time_min.copy(),  # tempo original; "time" pode ser deslocado pelo alinhamento de corridas
            "time": time_min,
            "tic": tic,
            "start": start,
            "count": count,
            "mass": mass,
            "inten": inten,
        },
    }


# ── Espectro médio ─────────────────────────────────────────────────────────

_MAX_BINS = 5_000_000


def _gather_points(arrays: dict, s0: int, s1: int):
    """Pontos (m/z, intensidade) de todos os scans de s0 a s1 (inclusive)."""
    start = arrays["start"][s0 : s1 + 1]
    count = arrays["count"][s0 : s1 + 1]
    a = int(start[0])
    b = int(start[-1] + count[-1])
    if len(start) == 1 or np.all(start[1:] == start[:-1] + count[:-1]):
        # scans contíguos nos vetores planos (caso normal): uma fatia só
        return arrays["mass"][a:b], arrays["inten"][a:b]
    idx = np.concatenate([np.arange(st, st + c) for st, c in zip(start, count)])
    return arrays["mass"][idx], arrays["inten"][idx]


def average_scans(arrays: dict, scan_start: int, scan_end: int, bin_width: float = 1.0, bg_start=None, bg_end=None) -> dict:
    """
    Espectro médio dos scans [scan_start, scan_end] com alinhamento de m/z por bins.

    Alinhamento: cada m/z é levado ao ponto mais próximo de uma grade de passo `bin_width`
    (1 = massa nominal, 0.1, 0.01 ...) e as intensidades que caem no mesmo bin são somadas
    dentro de cada scan. A média é sobre TODOS os scans do intervalo (um bin ausente em um scan
    conta como zero), de modo que o resultado seja a intensidade média por scan.

    Opcionalmente calcula a média de uma região de fundo na MESMA grade, devolvida em separado
    ("background") para o cliente poder ligar/desligar a subtração sem nova requisição.
    """
    if not (0.001 <= bin_width <= 10):
        raise ValueError("A largura do bin deve estar entre 0,001 e 10.")

    n_total = len(arrays["time"])

    def norm(a, b):
        a, b = sorted((int(a), int(b)))
        a, b = max(0, a), min(n_total - 1, b)
        if a > b:
            raise ValueError("Intervalo de scans fora do arquivo.")
        return a, b

    s0, s1 = norm(scan_start, scan_end)
    mass_a, int_a = _gather_points(arrays, s0, s1)
    n_a = s1 - s0 + 1

    has_bg = bg_start is not None and bg_end is not None
    if has_bg:
        b0, b1 = norm(bg_start, bg_end)
        mass_b, int_b = _gather_points(arrays, b0, b1)
        n_b = b1 - b0 + 1

    idx_a = np.rint(mass_a / bin_width).astype(np.int64)
    idx_b = np.rint(mass_b / bin_width).astype(np.int64) if has_bg else np.empty(0, dtype=np.int64)

    parts = [p for p in (idx_a, idx_b) if p.size]
    if not parts:
        return {
            "mz": [],
            "intensity": [],
            "background": None,
            "n_scans": n_a,
            "n_bg_scans": 0,
            "scan_start": s0,
            "scan_end": s1,
            "bg_start": None,
            "bg_end": None,
            "rt_start": float(arrays["time"][s0]),
            "rt_end": float(arrays["time"][s1]),
            "bin_width": bin_width,
        }
    lo = int(min(p.min() for p in parts))
    hi = int(max(p.max() for p in parts))
    size = hi - lo + 1
    if size > _MAX_BINS:
        raise ValueError("Largura de bin pequena demais para esta faixa de m/z.")

    avg = np.bincount(idx_a - lo, weights=int_a, minlength=size) / n_a
    if has_bg:
        bg = np.bincount(idx_b - lo, weights=int_b, minlength=size) / n_b
    else:
        bg = np.zeros(size)

    keep = np.nonzero((avg > 0) | (bg > 0))[0]
    return {
        "mz": np.round((keep + lo) * bin_width, 4).tolist(),
        "intensity": avg[keep].tolist(),
        "background": bg[keep].tolist() if has_bg else None,
        "n_scans": n_a,
        "n_bg_scans": n_b if has_bg else 0,
        "scan_start": s0,
        "scan_end": s1,
        "bg_start": b0 if has_bg else None,
        "bg_end": b1 if has_bg else None,
        "rt_start": float(arrays["time"][s0]),
        "rt_end": float(arrays["time"][s1]),
        "bin_width": bin_width,
    }


# ── Cromatograma de íon extraído (EIC) ─────────────────────────────────────


def extract_ion_chromatograms(arrays: dict, ions, tolerance: float = 0.5) -> list:
    """
    Para cada m/z em `ions`, soma por scan as intensidades dos pontos com |m/z - íon| <= tolerância.
    Devolve uma lista de séries, cada uma com um valor por scan (mesmo comprimento do TIC).
    """
    if not ions:
        raise ValueError("Informe ao menos um m/z.")
    if len(ions) > 16:
        raise ValueError("No máximo 16 íons por vez.")
    if not (0 < tolerance <= 5):
        raise ValueError("A tolerância deve estar entre 0 e 5 Da.")

    return [np.round(eic_array(arrays, ion, tolerance), 2).tolist() for ion in ions]


def eic_array(arrays: dict, ion: float, tolerance: float = 0.5):
    """EIC de um único m/z como array NumPy (um valor por scan)."""
    mass, inten = arrays["mass"], arrays["inten"]
    start = arrays["start"]
    end = start + arrays["count"]
    w = np.where(np.abs(mass - float(ion)) <= tolerance, inten, 0.0)
    cs = np.concatenate(([0.0], np.cumsum(w)))
    return np.maximum(cs[end] - cs[start], 0.0)  # a diferença de somas acumuladas pode dar -1e-4 por arredondamento
