from fastapi import FastAPI, UploadFile, File, HTTPException, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import numpy as np
import os
import re
import sys
import json
import subprocess
from datetime import datetime

from parsers.gcfid import read_gcfid_cdf
from parsers.gcms import read_gcms_cdf, average_scans, extract_ion_chromatograms, eic_array
from integration import integrate_peak, detect_peak
from spectral import parse_msp, compare_spectra, Library

app = FastAPI(title="Chromatography Viewer API")

# O app roda apenas em localhost. Aceitar qualquer origem ("*") permitiria que qualquer site aberto
# no navegador chamasse estas rotas (apagar/sobrescrever projetos, abrir diálogos de pasta).
_LOCAL_HOST_RE = re.compile(r"^(127\.0\.0\.1|localhost|\[::1\])(:\d+)?$")
_LOCAL_ORIGIN_RE = re.compile(r"^https?://(127\.0\.0\.1|localhost|\[::1\])(:\d+)?$")


@app.middleware("http")
async def only_local_requests(request: Request, call_next):
    host = request.headers.get("host", "")
    origin = request.headers.get("origin")
    # Host fora de localhost bloqueia DNS rebinding; Origin estranho bloqueia requisições de outros sites
    if not _LOCAL_HOST_RE.match(host) or (origin and not _LOCAL_ORIGIN_RE.match(origin)):
        return JSONResponse({"detail": "Requisição não permitida."}, status_code=403)
    response = await call_next(request)
    if not request.url.path.startswith("/api/"):
        # Arquivos do frontend: o navegador revalida a cada carga (evita ficar com app.js/style.css antigos)
        response.headers["Cache-Control"] = "no-cache"
    return response


app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=_LOCAL_ORIGIN_RE.pattern,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Diretório padrão para salvar projetos
DEFAULT_PROJECTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "projetos"))
os.makedirs(DEFAULT_PROJECTS_DIR, exist_ok=True)
PROJECTS_DIR = DEFAULT_PROJECTS_DIR

SAMPLES_CACHE: Dict[str, Dict[str, Any]] = {}

# ── GC-MS: dados completos dos scans (espectros são servidos sob demanda) ──
MS_CACHE: Dict[str, Dict[str, Any]] = {}
MS_SHIFT: Dict[str, float] = {}   # deslocamento de tempo (min) de cada arquivo GC-MS, vindo do alinhamento
# Cópia do .CDF original em disco, para os espectros sobreviverem a reinício do servidor / reabertura de projeto
GCMS_CACHE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "cache_gcms"))
os.makedirs(GCMS_CACHE_DIR, exist_ok=True)


def _ms_cache_path(sample_id: str) -> str:
    safe = re.sub(r"[^\w\-. ()]", "_", sample_id)
    return os.path.join(GCMS_CACHE_DIR, safe + ".cdf")


def _get_ms(sample_id: str):
    if sample_id in MS_CACHE:
        return MS_CACHE[sample_id]
    path = _ms_cache_path(sample_id)
    if os.path.exists(path):
        with open(path, "rb") as f:
            parsed = read_gcms_cdf(f.read(), sample_id + ".cdf")
        arrays = parsed["_arrays"]
        arrays["time"] = arrays["time0"] + MS_SHIFT.get(sample_id, 0.0)
        MS_CACHE[sample_id] = arrays
        return arrays
    return None


@app.post("/api/upload/gc-fid")
async def upload_gcfid_files(files: List[UploadFile] = File(...)):
    results = []
    errors = []

    for file in files:
        if not file.filename.lower().endswith(".cdf"):
            errors.append({"filename": file.filename, "error": "Extensão inválida (esperado .cdf)"})
            continue

        try:
            content = await file.read()
            parsed_data = read_gcfid_cdf(content, file.filename)

            sample_id = os.path.splitext(file.filename)[0].strip()
            parsed_data["sample_id"] = sample_id

            SAMPLES_CACHE[sample_id] = {
                "time": np.array(parsed_data["time"]),
                "intensity": np.array(parsed_data["intensity"]),
                "delay_time_sec": parsed_data["delay_time_sec"],
                "sampling_interval_sec": parsed_data["sampling_interval_sec"],
            }

            results.append(parsed_data)
        except Exception as e:
            errors.append({"filename": file.filename, "error": str(e)})

    return {"success": results, "errors": errors}


@app.post("/api/project/sync")
async def sync_project_samples(samples: List[Dict[str, Any]]):
    # O projeto aberto é a única fonte da verdade: substitui o cache em vez de só acrescentar
    SAMPLES_CACHE.clear()
    for s in samples:
        s_id = s["id"]
        SAMPLES_CACHE[s_id] = {
            "time": np.array(s["time"]),
            "intensity": np.array(s["intensity"]),
            "delay_time_sec": s.get("delay_time_sec", 0.0),
            "sampling_interval_sec": s.get("sampling_interval_sec", 0.02),
        }
    return {"status": "OK", "synced": len(samples)}


# ── ALINHAMENTO DE CORRIDAS (deslocamento do eixo do tempo) ──


class ShiftRequest(BaseModel):
    sample_id: str
    delta: float


@app.post("/api/samples/shift")
async def shift_fid_sample(req: ShiftRequest):
    """GC-FID: soma `delta` (min) ao tempo da amostra no cache (o cliente já deslocou o seu)."""
    d = SAMPLES_CACHE.get(req.sample_id)
    if d is None:
        raise HTTPException(status_code=404, detail="Amostra não encontrada no servidor.")
    d["time"] = d["time"] + req.delta
    return {"status": "OK"}


class MsShiftsRequest(BaseModel):
    shifts: Dict[str, float]


@app.post("/api/gcms/shifts")
async def set_gcms_shifts(req: MsShiftsRequest):
    """GC-MS: define o deslocamento ABSOLUTO (min) de cada arquivo; vale também após reiniciar o servidor."""
    for sid, shift in req.shifts.items():
        MS_SHIFT[sid] = float(shift)
        if sid in MS_CACHE:
            MS_CACHE[sid]["time"] = MS_CACHE[sid]["time0"] + float(shift)
    return {"status": "OK"}


class AverageSpectrumRequest(BaseModel):
    sample_id: str
    scan_start: int
    scan_end: int
    bin_width: float = 1.0
    bg_start: Optional[int] = None
    bg_end: Optional[int] = None


@app.post("/api/gcms/average")
async def get_gcms_average(req: AverageSpectrumRequest):
    data = _get_ms(req.sample_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Dados do GC-MS não encontrados no servidor. Reimporte o arquivo .CDF.")
    try:
        return average_scans(data, req.scan_start, req.scan_end, req.bin_width, req.bg_start, req.bg_end)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


class EicRequest(BaseModel):
    sample_id: str
    ions: List[float]
    tol: float = 0.5


@app.post("/api/gcms/eic")
async def get_gcms_eic(req: EicRequest):
    data = _get_ms(req.sample_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Dados do GC-MS não encontrados no servidor. Reimporte o arquivo .CDF.")
    try:
        series = extract_ion_chromatograms(data, req.ions, req.tol)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    return {"ions": req.ions, "tol": req.tol, "series": series}


@app.post("/api/gcms/release")
async def release_gcms_memory():
    # Libera a memória dos espectros GC-MS; as cópias em cache_gcms/ continuam e são relidas sob demanda
    MS_CACHE.clear()
    return {"status": "OK"}


# ── GC-MS: BIBLIOTECAS MSP E COMPARAÇÃO DE ESPECTROS ──
LIBS: Dict[str, Library] = {}
LIB_DIR = os.path.join(GCMS_CACHE_DIR, "libs")
os.makedirs(LIB_DIR, exist_ok=True)


def _decode_text(content: bytes) -> str:
    try:
        return content.decode("utf-8")
    except UnicodeDecodeError:
        return content.decode("latin-1")


def _lib_id(filename: str) -> str:
    return re.sub(r"[^\w\-. ()]", "_", os.path.splitext(os.path.basename(filename))[0].strip()) or "biblioteca"


def _load_libs_from_disk():
    """Bibliotecas enviadas antes (mesmo após reiniciar o servidor) ficam em cache_gcms/libs/."""
    for f in os.listdir(LIB_DIR):
        lid = os.path.splitext(f)[0]
        if f.lower().endswith(".msp") and lid not in LIBS:
            with open(os.path.join(LIB_DIR, f), "rb") as fh:
                entries = parse_msp(_decode_text(fh.read()))
            if entries:
                LIBS[lid] = Library(lid, entries)


@app.post("/api/gcms/library/upload")
async def upload_gcms_library(files: List[UploadFile] = File(...)):
    results, errors = [], []
    for file in files:
        try:
            content = await file.read()
            entries = parse_msp(_decode_text(content))
            if not entries:
                raise ValueError("Nenhum espectro encontrado (formato MSP esperado: Name:, Num Peaks:, pares m/z intensidade).")
            lid = _lib_id(file.filename)
            LIBS[lid] = Library(lid, entries)
            with open(os.path.join(LIB_DIR, lid + ".msp"), "wb") as f:
                f.write(content)
            results.append({"lib_id": lid, "name": lid, "n_entries": len(entries)})
        except Exception as e:
            errors.append({"filename": file.filename, "error": str(e)})
    return {"success": results, "errors": errors}


@app.get("/api/gcms/library/list")
async def list_gcms_libraries():
    _load_libs_from_disk()
    return {"libraries": [{"lib_id": k, "name": v.name, "n_entries": len(v)} for k, v in LIBS.items()]}


class LibSearchRequest(BaseModel):
    lib_id: Optional[str] = None   # vazio = todas as bibliotecas
    mz: List[float]
    intensity: List[float]
    top: int = 10


@app.post("/api/gcms/library/search")
async def search_gcms_library(req: LibSearchRequest):
    _load_libs_from_disk()
    if req.lib_id and req.lib_id not in LIBS:
        raise HTTPException(status_code=404, detail="Biblioteca não encontrada.")
    libs = [(req.lib_id, LIBS[req.lib_id])] if req.lib_id else list(LIBS.items())
    hits = []
    for lid, lib in libs:
        for h in lib.search(req.mz, req.intensity, max(1, min(req.top, 50))):
            h["lib_id"] = lid
            hits.append(h)
    hits.sort(key=lambda h: -h["mf"])
    return {"hits": hits[: max(1, min(req.top, 50))]}


class LibFindRequest(BaseModel):
    lib_id: Optional[str] = None
    q: str


@app.post("/api/gcms/library/find")
async def find_gcms_library(req: LibFindRequest):
    _load_libs_from_disk()
    libs = [(req.lib_id, LIBS[req.lib_id])] if req.lib_id and req.lib_id in LIBS else list(LIBS.items())
    out = []
    for lid, lib in libs:
        for h in lib.find(req.q):
            h["lib_id"] = lid
            out.append(h)
    return {"hits": out[:30]}


@app.get("/api/gcms/library/entry/{lib_id}/{index}")
async def get_gcms_library_entry(lib_id: str, index: int):
    _load_libs_from_disk()
    lib = LIBS.get(lib_id)
    if lib is None or index < 0 or index >= len(lib):
        raise HTTPException(status_code=404, detail="Entrada não encontrada.")
    return lib.entry(index)


class CompareRequest(BaseModel):
    mz: List[float]
    intensity: List[float]
    ref_mz: List[float]
    ref_intensity: List[float]


@app.post("/api/gcms/compare")
async def compare_gcms_spectra(req: CompareRequest):
    return compare_spectra(req.mz, req.intensity, req.ref_mz, req.ref_intensity)


# ── GC-MS ──


@app.post("/api/upload/gc-ms")
async def upload_gcms_files(files: List[UploadFile] = File(...)):
    results = []
    errors = []

    for file in files:
        if not file.filename.lower().endswith(".cdf"):
            errors.append({"filename": file.filename, "error": "Extensão inválida (esperado .cdf)"})
            continue

        try:
            content = await file.read()
            parsed = read_gcms_cdf(content, file.filename)

            sample_id = os.path.splitext(file.filename)[0].strip()
            MS_CACHE[sample_id] = parsed.pop("_arrays")
            MS_SHIFT[sample_id] = 0.0   # arquivo novo/reimportado volta ao tempo original
            with open(_ms_cache_path(sample_id), "wb") as f:
                f.write(content)

            parsed["sample_id"] = sample_id
            results.append(parsed)
        except Exception as e:
            errors.append({"filename": file.filename, "error": str(e)})

    return {"success": results, "errors": errors}


@app.get("/api/gcms/spectrum/{sample_id}/{scan}")
async def get_gcms_spectrum(sample_id: str, scan: int):
    data = _get_ms(sample_id)
    if data is None:
        raise HTTPException(status_code=404, detail="Dados do GC-MS não encontrados no servidor. Reimporte o arquivo .CDF.")

    n_scans = len(data["time"])
    if scan < 0 or scan >= n_scans:
        raise HTTPException(status_code=400, detail=f"Scan fora do intervalo (0 a {n_scans - 1}).")

    a = int(data["start"][scan])
    b = a + int(data["count"][scan])
    return {
        "scan": scan,
        "rt_min": float(data["time"][scan]),
        "tic": float(data["tic"][scan]),
        "mz": np.round(data["mass"][a:b], 4).tolist(),
        "intensity": data["inten"][a:b].tolist(),
    }


# ── ROTAS DE GERENCIAMENTO DO DIRETÓRIO DE PROJETOS ──


@app.get("/api/projects/info")
async def get_projects_info():
    global PROJECTS_DIR
    entries = []
    if os.path.exists(PROJECTS_DIR):
        for f in os.listdir(PROJECTS_DIR):
            if f.lower().endswith(".json"):
                f_path = os.path.join(PROJECTS_DIR, f)
                stat = os.stat(f_path)
                entries.append(
                    (
                        stat.st_mtime,
                        {
                            "filename": f,
                            "name": os.path.splitext(f)[0],
                            "size_kb": round(stat.st_size / 1024, 1),
                            "modified": datetime.fromtimestamp(stat.st_mtime).strftime("%d/%m/%Y %H:%M"),
                        },
                    )
                )
    entries.sort(key=lambda e: e[0], reverse=True)  # mais recente primeiro, pela data numérica
    return {"directory": PROJECTS_DIR, "projects": [e[1] for e in entries]}


def _project_file_path(filename: str) -> str:
    """Resolve um nome de projeto dentro de PROJECTS_DIR, recusando caminhos (..\\, /, etc.)."""
    if (
        not filename
        or filename != os.path.basename(filename.replace("\\", "/"))
        or not filename.lower().endswith(".json")
    ):
        raise HTTPException(status_code=400, detail="Nome de arquivo inválido.")
    base = os.path.normcase(os.path.realpath(PROJECTS_DIR))
    path = os.path.realpath(os.path.join(PROJECTS_DIR, filename))
    if os.path.normcase(os.path.dirname(path)) != base:
        raise HTTPException(status_code=400, detail="Nome de arquivo inválido.")
    return path


# Diálogos nativos (Explorer) rodam em um processo separado: o tkinter exige a thread principal,
# e as rotas síncronas do FastAPI rodam em threads de trabalho (podia travar ou falhar em silêncio).
_TK_DIALOG_SCRIPT = r"""
import sys, json
import tkinter as tk
from tkinter import filedialog

kind, initial_dir, initial_file = sys.argv[1], sys.argv[2], sys.argv[3]
root = tk.Tk()
# Janela-dona 1x1 e transparente (em vez de withdraw): um processo em segundo plano não consegue
# trazer um diálogo para a frente sozinho, mas um diálogo de dono "topmost" aparece sobre o navegador.
try:
    root.attributes("-topmost", True)
    root.attributes("-alpha", 0.0)
except Exception:
    pass
root.geometry("1x1+0+0")
root.update()
root.lift()
root.focus_force()
if kind == "dir":
    sel = filedialog.askdirectory(parent=root, initialdir=initial_dir, title="Selecione a Pasta de Projetos")
else:
    sel = filedialog.asksaveasfilename(
        parent=root,
        initialdir=initial_dir,
        initialfile=initial_file,
        defaultextension=".json",
        filetypes=[("Arquivo de Projeto GC (*.json)", "*.json"), ("Todos os Arquivos (*.*)", "*.*")],
        title="Salvar Projeto Como",
    )
root.destroy()
print(json.dumps({"path": sel if isinstance(sel, str) else ""}))
"""


def _run_tk_dialog(kind: str, initial_dir: str, initial_file: str = ""):
    """Retorna (status, caminho, mensagem): status OK | CANCELLED | ERROR (erro real, não confundir com cancelar)."""
    try:
        flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        proc = subprocess.run(
            [sys.executable, "-c", _TK_DIALOG_SCRIPT, kind, initial_dir or "", initial_file or ""],
            capture_output=True,
            text=True,
            timeout=900,
            creationflags=flags,
        )
        if proc.returncode != 0:
            msg = (proc.stderr or "").strip().splitlines()
            return "ERROR", "", msg[-1] if msg else "Falha ao abrir a janela de seleção."
        path = json.loads(proc.stdout.strip().splitlines()[-1]).get("path", "").strip()
        return ("OK", path, "") if path else ("CANCELLED", "", "")
    except Exception as e:
        return "ERROR", "", str(e)


@app.post("/api/projects/pick-directory")
def pick_directory_native():
    global PROJECTS_DIR
    status, path, message = _run_tk_dialog("dir", PROJECTS_DIR)
    if status == "OK":
        PROJECTS_DIR = os.path.normpath(path)
        return {"status": "OK", "directory": PROJECTS_DIR}
    if status == "ERROR":
        return {"status": "ERROR", "message": message, "directory": PROJECTS_DIR}
    return {"status": "CANCELLED", "directory": PROJECTS_DIR}


class PickSavePathRequest(BaseModel):
    default_name: Optional[str] = None
    initial_dir: Optional[str] = None


# Janela nativa "Salvar Como" (local e nome do arquivo)
@app.post("/api/projects/pick-save-path")
def pick_save_path_native(req: PickSavePathRequest):
    global PROJECTS_DIR
    init_dir = req.initial_dir or PROJECTS_DIR
    if not init_dir or not os.path.exists(init_dir):
        init_dir = os.path.expanduser("~")

    gen_name = req.default_name.strip() if req.default_name else f"Projeto_{datetime.now().strftime('%Y-%m-%d')}"
    if not gen_name.lower().endswith(".json"):
        gen_name += ".json"

    status, selected_file, message = _run_tk_dialog("save", init_dir, gen_name)
    if status == "OK":
        norm_path = os.path.normpath(selected_file)
        new_dir = os.path.dirname(norm_path)
        if new_dir:
            PROJECTS_DIR = new_dir
        return {
            "status": "OK",
            "filepath": norm_path,
            "directory": new_dir,
            "filename": os.path.basename(norm_path),
        }
    if status == "ERROR":
        return {"status": "ERROR", "message": message}
    return {"status": "CANCELLED"}


class SetDirectoryRequest(BaseModel):
    directory: str


@app.post("/api/projects/set-directory")
async def set_projects_directory(req: SetDirectoryRequest):
    global PROJECTS_DIR
    path = os.path.normpath(req.directory.strip())
    os.makedirs(path, exist_ok=True)
    PROJECTS_DIR = path
    return {"status": "OK", "directory": PROJECTS_DIR}


class SaveProjectRequest(BaseModel):
    name: Optional[str] = None
    filepath: Optional[str] = None
    project_data: Dict[str, Any]


@app.post("/api/projects/save")
async def save_project_to_dir(req: SaveProjectRequest):
    global PROJECTS_DIR

    if req.filepath and req.filepath.strip():
        file_path = os.path.normpath(req.filepath.strip())
        if not file_path.lower().endswith(".json"):
            file_path += ".json"
        target_dir = os.path.dirname(file_path)
        if target_dir:
            os.makedirs(target_dir, exist_ok=True)
            PROJECTS_DIR = target_dir
        safe_name = os.path.basename(file_path)
    else:
        raw_name = req.name.strip() if req.name and req.name.strip() else f"Projeto_{datetime.now().strftime('%Y-%m-%d')}"
        safe_name = raw_name.replace(" ", "_")
        if not safe_name.lower().endswith(".json"):
            safe_name += ".json"
        os.makedirs(PROJECTS_DIR, exist_ok=True)
        file_path = os.path.join(PROJECTS_DIR, safe_name)

    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(req.project_data, f, indent=2, ensure_ascii=False)

    return {"status": "OK", "filename": safe_name, "path": file_path, "directory": PROJECTS_DIR}


@app.get("/api/projects/load/{filename}")
async def load_project_from_dir(filename: str):
    file_path = _project_file_path(filename)
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="Projeto não encontrado.")

    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    return {"status": "OK", "data": data, "filename": filename, "filepath": file_path}


@app.delete("/api/projects/delete/{filename}")
async def delete_project_from_dir(filename: str):
    file_path = _project_file_path(filename)
    if os.path.exists(file_path):
        os.remove(file_path)
    return {"status": "OK"}


# ── SINAIS INTEGRÁVEIS ──
# Um "sinal" é uma amostra de GC-FID (id = nome do arquivo) ou uma curva do GC-MS no formato
# "msig|<arquivo>|<íon ou TIC>|<tolerância>": o EIC de um m/z (ou o TIC) de um arquivo GC-MS.


def parse_signal_id(sample_id: str):
    """Retorna (ms_id, ion|None, tol) para ids de GC-MS; None para amostras de FID."""
    if not sample_id.startswith("msig|"):
        return None
    try:
        ms_id, ion, tol = sample_id[5:].rsplit("|", 2)
        return ms_id, (None if ion == "TIC" else float(ion)), float(tol)
    except ValueError:
        return None


def get_signal_arrays(sample_id: str):
    """(tempo em min, intensidade) do sinal, ou None se os dados não estiverem no servidor."""
    ms = parse_signal_id(sample_id)
    if ms is not None:
        ms_id, ion, tol = ms
        data = _get_ms(ms_id)
        if data is None:
            return None
        if ion is None:
            return data["time"], data["tic"]
        return data["time"], eic_array(data, ion, tol)
    d = SAMPLES_CACHE.get(sample_id)
    return (d["time"], d["intensity"]) if d else None


def signal_info(sample_id: str):
    ms = parse_signal_id(sample_id)
    if ms is None:
        return None
    return {"ms_id": ms[0], "ion": ms[1], "tol": ms[2]}


class ForcedIntegrationRequest(BaseModel):
    sample_ids: List[str]
    peak_name: str
    baseline_type: str
    start_idx: int
    end_idx: int
    gl_start_idx: Optional[int] = None
    gl_end_idx: Optional[int] = None
    anchor_idx: Optional[int] = None


@app.post("/api/integrate/forced")
async def calculate_forced_integration(req: ForcedIntegrationRequest):
    integrated_results = []

    for s_id in req.sample_ids:
        arrays = get_signal_arrays(s_id)
        if arrays is None:
            continue
        res = integrate_peak(arrays[0], arrays[1], req.baseline_type, req.start_idx, req.end_idx,
                             req.gl_start_idx, req.gl_end_idx, req.anchor_idx)
        if res is None:
            continue
        res["sample_id"] = s_id
        res["peak_name"] = req.peak_name
        res["signal"] = signal_info(s_id)
        integrated_results.append(res)

    return {"status": "OK", "peaks": integrated_results}


# ── MÉTODO DE COMPOSTOS: detecção automática do pico na janela de cada composto ──


class DetectItem(BaseModel):
    signal_id: str
    compound_id: str
    rt: float
    win: float
    baseline_type: str = "Valley-to-Valley"


class DetectRequest(BaseModel):
    items: List[DetectItem]
    snr_min: float = 3.0


@app.post("/api/compounds/detect")
async def detect_compounds(req: DetectRequest):
    results = []
    cache: Dict[str, Any] = {}
    for it in req.items:
        if it.signal_id not in cache:
            cache[it.signal_id] = get_signal_arrays(it.signal_id)
        arrays = cache[it.signal_id]
        if arrays is None:
            results.append({"compound_id": it.compound_id, "signal_id": it.signal_id, "found": False,
                            "reason": "dados da corrida não estão no servidor"})
            continue
        peak, info = detect_peak(arrays[0], arrays[1], it.rt, max(it.win, 0.001), it.baseline_type, req.snr_min)
        if peak is None:
            results.append({"compound_id": it.compound_id, "signal_id": it.signal_id, "found": False, "reason": info})
            continue
        peak["sample_id"] = it.signal_id
        peak["signal"] = signal_info(it.signal_id)
        results.append({"compound_id": it.compound_id, "signal_id": it.signal_id, "found": True,
                        "snr": info["snr"], "peak": peak})
    return {"status": "OK", "results": results}


frontend_dir = os.path.join(os.path.dirname(__file__), "..", "frontend")
app.mount("/", StaticFiles(directory=frontend_dir, html=True), name="frontend")

if __name__ == "__main__":
    import uvicorn

    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)