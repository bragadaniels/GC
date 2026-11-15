"""Comparação de espectros de massas: leitura de MSP, nota de semelhança (estilo NIST) e busca em biblioteca."""

import re
import numpy as np
from scipy import sparse

_NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
MZ_POWER = 3.0  # peso do m/z na nota (íons pesados contam mais) — convenção NIST
INT_POWER = 0.6  # peso da intensidade (reduz a dominância do pico base)
MAX_MZ = 2000


# ── Leitura de MSP ─────────────────────────────────────────────────────────


def parse_msp(text: str) -> list:
    """
    Lê um arquivo MSP. Cada entrada: campos "Chave: valor", "Num Peaks: N" e depois os pares m/z intensidade
    (separados por espaço, tab, vírgula ou ';', vários por linha ou um por linha), até uma linha em branco.
    Retorna [{name, meta, mz, intensity}].
    """
    entries = []
    cur = None
    in_peaks = False
    peak_text = []

    def flush():
        nonlocal cur, in_peaks, peak_text
        if cur is not None:
            nums = [float(x) for x in _NUM.findall(" ".join(peak_text))]
            if len(nums) >= 2:
                pairs = np.array(nums[: len(nums) // 2 * 2]).reshape(-1, 2)
                cur["mz"], cur["intensity"] = pairs[:, 0], pairs[:, 1]
                if cur.get("name"):
                    entries.append(cur)
        cur, in_peaks, peak_text = None, False, []

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            flush()
            continue
        if in_peaks and (line[0].isdigit() or line[0] in "+-."):
            peak_text.append(line)
            continue
        m = re.match(r"^([A-Za-z_ ]+?)\s*:\s*(.*)$", line)
        if m:
            key, val = m.group(1).strip().lower().replace(" ", "_"), m.group(2).strip()
            if key == "name":
                flush()
                cur = {"name": val, "meta": {}}
            elif cur is None:
                continue
            elif key in ("num_peaks", "numpeaks"):
                in_peaks = True
                rest = val  # alguns arquivos já trazem os picos na mesma linha depois do número
                tail = _NUM.findall(rest)
                if len(tail) > 1:
                    peak_text.append(" ".join(tail[1:]))
            else:
                cur["meta"][key] = val
            continue
        if cur is not None and in_peaks:
            peak_text.append(line)
    flush()
    return entries


# ── Vetores e nota ─────────────────────────────────────────────────────────


def to_nominal(mz, inten):
    """Agrupa em massa nominal (arredonda o m/z) somando as intensidades. Retorna (índices únicos, soma)."""
    mz = np.asarray(mz, dtype=np.float64)
    inten = np.asarray(inten, dtype=np.float64)
    idx = np.rint(mz).astype(np.int64)
    ok = (idx >= 1) & (idx <= MAX_MZ) & (inten > 0)
    idx, inten = idx[ok], inten[ok]
    if idx.size == 0:
        return idx, inten
    uniq, inv = np.unique(idx, return_inverse=True)
    return uniq, np.bincount(inv, weights=inten)


def _weights(idx, vals):
    return (idx.astype(np.float64) ** MZ_POWER) * (vals**INT_POWER)


def _dense(idx, w, size):
    v = np.zeros(size)
    v[idx] = w
    return v


def compare_spectra(q_mz, q_int, r_mz, r_int) -> dict:
    """
    Nota de semelhança 0–999 entre uma consulta e uma referência (cosseno com pesos m/z^3 · I^0,6):
      mf  = direta (todos os íons);  rmf = reversa (só os íons presentes na referência).
    Também devolve os vetores alinhados em massa nominal (normalizados a 100) para o gráfico espelhado.
    """
    qi, qv = to_nominal(q_mz, q_int)
    ri, rv = to_nominal(r_mz, r_int)
    if qi.size == 0 or ri.size == 0:
        return {
            "mf": 0.0,
            "rmf": 0.0,
            "n_common": 0,
            "n_query": int(qi.size),
            "n_ref": int(ri.size),
            "mz": [],
            "query": [],
            "ref": [],
        }
    size = int(max(qi.max(), ri.max())) + 1
    wq, wr = _dense(qi, _weights(qi, qv), size), _dense(ri, _weights(ri, rv), size)
    nq, nr = np.linalg.norm(wq), np.linalg.norm(wr)
    mf = 999.0 * float(wq @ wr) / (nq * nr) if nq > 0 and nr > 0 else 0.0

    mask = wr > 0
    wqm = wq * mask
    nqm = np.linalg.norm(wqm)
    rmf = 999.0 * float(wqm @ wr) / (nqm * nr) if nqm > 0 and nr > 0 else 0.0

    union = np.union1d(qi, ri)
    q_plot = _dense(qi, qv, size)[union]
    r_plot = _dense(ri, rv, size)[union]
    return {
        "mf": min(mf, 999.0),
        "rmf": min(rmf, 999.0),
        "n_common": int(np.sum((wq > 0) & (wr > 0))),
        "n_query": int(qi.size),
        "n_ref": int(ri.size),
        "mz": union.astype(float).tolist(),
        "query": (q_plot / q_plot.max() * 100.0).tolist(),
        "ref": (r_plot / r_plot.max() * 100.0).tolist(),
    }


# ── Biblioteca ─────────────────────────────────────────────────────────────


class Library:
    def __init__(self, name: str, entries: list):
        self.name = name
        self.names = [e["name"] for e in entries]
        self.meta = [e["meta"] for e in entries]
        self.size = MAX_MZ + 1
        rows, cols, vals = [], [], []
        for r, e in enumerate(entries):
            idx, v = to_nominal(e["mz"], e["intensity"])
            rows.append(np.full(idx.size, r))
            cols.append(idx)
            vals.append(v)
        n = len(entries)
        if rows:
            rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
        else:
            rows = cols = np.empty(0, dtype=np.int64)
            vals = np.empty(0)
        # intensidades por massa nominal (para recuperar o espectro de uma entrada)
        self.raw = sparse.csr_matrix((vals.astype(np.float32), (rows, cols)), shape=(n, self.size))
        # versão ponderada e normalizada (norma 1) para o cosseno; máscara dos íons de cada entrada
        w = _weights(cols, vals)
        W = sparse.csr_matrix((w, (rows, cols)), shape=(n, self.size))
        norms = np.sqrt(np.asarray(W.multiply(W).sum(axis=1)).ravel())
        norms[norms == 0] = 1.0
        self.W = sparse.diags(1.0 / norms) @ W
        self.mask = (W != 0).astype(np.float64).tocsr()

    def __len__(self):
        return len(self.names)

    def entry(self, i: int) -> dict:
        row = self.raw.getrow(i)
        return {
            "name": self.names[i],
            "meta": self.meta[i],
            "mz": row.indices.astype(float).tolist(),
            "intensity": row.data.astype(float).tolist(),
        }

    def search(self, q_mz, q_int, top: int = 10) -> list:
        qi, qv = to_nominal(q_mz, q_int)
        qi_ok = qi < self.size
        qi, qv = qi[qi_ok], qv[qi_ok]
        if qi.size == 0 or len(self) == 0:
            return []
        wq = _dense(qi, _weights(qi, qv), self.size)
        nq = np.linalg.norm(wq)
        num = self.W @ wq
        fwd = num / nq
        denom = np.sqrt(self.mask @ (wq**2))
        rev = np.divide(num, denom, out=np.zeros_like(num), where=denom > 0)
        common = (self.mask @ (wq > 0).astype(np.float64)).astype(int)
        order = np.argsort(-fwd)[:top]
        return [
            {
                "index": int(i),
                "name": self.names[i],
                "meta": self.meta[i],
                "mf": min(999.0 * float(fwd[i]), 999.0),
                "rmf": min(999.0 * float(rev[i]), 999.0),
                "n_common": int(common[i]),
            }
            for i in order
        ]

    def find(self, text: str, limit: int = 30) -> list:
        t = text.strip().lower()
        out = []
        for i, n in enumerate(self.names):
            if t in n.lower():
                out.append({"index": i, "name": n, "meta": self.meta[i]})
                if len(out) >= limit:
                    break
        return out
