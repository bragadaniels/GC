"""Matemática de integração de picos (compartilhada por FID, EIC/TIC do GC-MS e pelo método de compostos)."""

import numpy as np
from scipy.integrate import trapezoid


def integrate_peak(time, intensity, baseline_type, start_idx, end_idx, gl_start_idx=None, gl_end_idx=None, anchor_idx=None):
    """
    Integra um pico entre start_idx e end_idx (índices de pontos) com a linha de base escolhida.
    `time` em minutos; a área sai em (unidade do sinal · s), com precisão total (sem arredondar).
    Retorna None se o intervalo for inválido.
    """
    time = np.asarray(time, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)
    n_points = len(time)

    start = max(0, min(int(start_idx), n_points - 1))
    end = max(0, min(int(end_idx), n_points - 1))
    if start >= end:
        return None

    x_peak = time[start : end + 1]
    y_peak = intensity[start : end + 1]

    if baseline_type == "Valley-to-Valley":
        y_base = np.linspace(y_peak[0], y_peak[-1], len(x_peak))
        bl_x = [float(x_peak[0]), float(x_peak[-1])]
        bl_y = [float(y_peak[0]), float(y_peak[-1])]

    elif baseline_type == "Drop Line":
        gl_s = max(0, min(gl_start_idx if gl_start_idx is not None else start, n_points - 1))
        gl_e = max(0, min(gl_end_idx if gl_end_idx is not None else end, n_points - 1))
        if gl_s >= gl_e:
            gl_s, gl_e = start, end

        x_cluster = time[gl_s : gl_e + 1]
        y_cluster_base = np.linspace(intensity[gl_s], intensity[gl_e], len(x_cluster))
        y_base = np.interp(x_peak, x_cluster, y_cluster_base)
        bl_x = [float(time[gl_s]), float(time[gl_e])]
        bl_y = [float(intensity[gl_s]), float(intensity[gl_e])]

    elif baseline_type == "Tangent Skim":
        anc = max(0, min(anchor_idx if anchor_idx is not None else end, n_points - 1))
        x_anc = time[anc]
        y_anc = intensity[anc]
        x_span = [time[start], x_anc]
        y_span = [intensity[start], y_anc]
        y_base = np.interp(x_peak, x_span, y_span)
        bl_x = [float(time[start]), float(x_anc)]
        bl_y = [float(intensity[start]), float(y_anc)]
    else:
        y_base = np.linspace(y_peak[0], y_peak[-1], len(x_peak))
        bl_x = [float(x_peak[0]), float(x_peak[-1])]
        bl_y = [float(y_peak[0]), float(y_peak[-1])]

    y_above = np.maximum(y_peak - y_base, 0.0)

    # Integral com tempo real em segundos e precisão TOTAL (sem round)
    time_seconds = time * 60.0
    area_calc = float(trapezoid(y_above, x=time_seconds[start : end + 1]))

    apex_local_idx = int(np.argmax(y_above))
    apex_global_idx = start + apex_local_idx
    apex_rt_min = float(time[apex_global_idx])
    apex_intensity = float(intensity[apex_global_idx])
    height_calc = float(y_above[apex_local_idx])

    poly_x = np.concatenate([x_peak, x_peak[::-1]]).tolist()
    poly_y = np.concatenate([y_peak, y_base[::-1]]).tolist()

    return {
        "baseline_type": baseline_type,
        "area": area_calc,
        "height": height_calc,
        "rt_min": apex_rt_min,
        "apex": {"x": apex_rt_min, "y": apex_intensity},
        "baseline_line": {"x": bl_x, "y": bl_y},
        "polygon": {"x": poly_x, "y": poly_y},
    }


# ── Detecção automática (método de compostos) ──────────────────────────────


def estimate_noise(y) -> float:
    """Ruído robusto (desvio) pela mediana das diferenças entre pontos vizinhos; ignora os picos."""
    y = np.asarray(y, dtype=np.float64)
    if len(y) < 5:
        return 0.0
    d = np.diff(y)
    mad = np.median(np.abs(d - np.median(d)))
    return float(1.4826 * mad / np.sqrt(2.0))


def detect_peak(time, y, rt, win, baseline_type="Valley-to-Valley", snr_min=3.0):
    """
    Procura um pico no tR esperado (rt ± win, em minutos) e o integra.

    - ápice = máximo do sinal na janela (tem de ficar dentro dela, não na borda);
    - S/N = (ápice − base local) / ruído do cromatograma inteiro, mínimo `snr_min`;
    - limites: a partir do ápice, desce até o vale (o sinal volta a subir) ou até cair a ~1% da altura,
      podendo passar da janela até 2× a largura dela (para picos largos).
    Retorna (peak_dict, info) ou (None, motivo).
    """
    time = np.asarray(time, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = len(time)
    lo = int(np.searchsorted(time, rt - win, side="left"))
    hi = int(np.searchsorted(time, rt + win, side="right"))  # exclusivo
    if hi - lo < 3:
        return None, "janela fora do cromatograma"

    noise = estimate_noise(y)
    # suavização leve: o ápice e o caminho até o vale usam a curva suavizada (um pico de ruído não vira "ápice")
    k = 5
    ys = np.convolve(y, np.ones(k) / k, mode="same")
    noise_s = noise / np.sqrt(k)

    apex = lo + int(np.argmax(ys[lo:hi]))
    if apex == lo or apex == hi - 1:
        return None, "sem pico dentro da janela (máximo na borda)"

    ext_lo = max(0, int(np.searchsorted(time, rt - 2 * win, side="left")))
    ext_hi = min(n - 1, int(np.searchsorted(time, rt + 2 * win, side="right")) - 1)

    # base local: mínimos de cada lado na faixa ESTENDIDA (a janela pode estar toda em cima de um pico largo)
    base_guess = 0.5 * (ys[ext_lo : apex + 1].min() + ys[apex : ext_hi + 1].min())
    height = ys[apex] - base_guess
    snr = height / max(noise, 1e-12 * max(1.0, float(np.max(np.abs(y)))))
    if snr < snr_min:
        return None, f"S/N baixo ({snr:.1f} < {snr_min:g})"

    floor = base_guess + 0.01 * height
    # O vale é o mínimo corrente, confirmado quando o sinal volta a subir de forma significativa
    # (robusto ao ruído proporcional ao sinal, comum em MS, que faz a flanco "serrilhar").
    # ruído medido no corpo do pico (resíduo contra a curva suavizada): captura ruído proporcional ao sinal
    resid = (y - ys)[ext_lo : ext_hi + 1]
    ys_seg = ys[ext_lo : ext_hi + 1]
    body = ys_seg > floor
    rel_noise = float(np.std(resid[body] / ys_seg[body])) if body.sum() >= 5 else 0.0
    rel_noise = min(rel_noise, 0.5)
    rise_abs = max(5 * noise_s, 0.03 * height)
    rise_rel = 5 * rel_noise / np.sqrt(k)  # fração do nível local do sinal
    valley_zone = base_guess + 0.85 * height  # só se procura vale abaixo de 85% da altura do pico

    def walk(step, limit):
        i = apex
        m_idx, m_val = apex, ys[apex]
        while True:
            nxt = i + step
            if (step < 0 and nxt < limit) or (step > 0 and nxt > limit):
                break
            i = nxt
            if ys[i] < m_val:
                m_idx, m_val = i, ys[i]
            if ys[i] <= floor:
                m_idx = i
                break
            if ys[i] < valley_zone and ys[i] - m_val > max(rise_abs, rise_rel * ys[i]):
                break
        return m_idx

    start = walk(-1, ext_lo)
    end = walk(+1, ext_hi)
    if end - start < 2:
        return None, "pico estreito demais"

    peak = integrate_peak(time, y, baseline_type, start, end)
    if peak is None:
        return None, "intervalo inválido"
    return peak, {"snr": float(snr), "start_idx": int(start), "end_idx": int(end)}
