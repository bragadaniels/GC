// ── SELEÇÃO DE ELEMENTOS DO DOM ──
const sidebar = document.getElementById('sidebar');
const resizer = document.getElementById('resizer');
const collapseBtn = document.getElementById('collapse-btn');
const expandBtn = document.getElementById('expand-btn');
const fileInput = document.getElementById('file-input');
const uploadBtn = document.getElementById('upload-btn');
const clearBtn = document.getElementById('clear-btn');
const themeBtn = document.getElementById('theme-btn');
const themeIcon = document.getElementById('theme-icon');
const sampleList = document.getElementById('sample-list');
const sampleCount = document.getElementById('sample-count');
const plotDiv = document.getElementById('plot');
const chartArea = document.getElementById('chart-area');

// Painel Direito de Projetos
const rightSidebar = document.getElementById('right-sidebar');
const rightResizer = document.getElementById('right-resizer');
const toggleProjectsPanelBtn = document.getElementById('toggle-projects-panel-btn');
const closeProjectsBtn = document.getElementById('close-projects-btn');
const currentDirDisplay = document.getElementById('current-dir-display');
const currentProjectDisplay = document.getElementById('current-project-display');
const btnChangeDirectory = document.getElementById('btn-change-directory');
const btnSaveProject = document.getElementById('btn-save-project');
const btnSaveAsProject = document.getElementById('btn-save-as-project');
const btnLoadExternalProject = document.getElementById('btn-load-external-project');
const btnCloseProject = document.getElementById('btn-close-project');
const externalProjectInput = document.getElementById('external-project-input');
const projectsList = document.getElementById('projects-list');
const projectsCount = document.getElementById('projects-count');

// Integração & Banner
const integrateBtn = document.getElementById('integrate-btn');
const toggleTableBtn = document.getElementById('toggle-table-btn');
const notesBtn = document.getElementById('notes-btn');
const tempRampBtn = document.getElementById('temp-ramp-btn');
const integrationBanner = document.getElementById('integration-banner');
const bannerText = document.getElementById('banner-text');
const bannerCancelBtn = document.getElementById('banner-cancel-btn');

// Bloco de Notas
const notesModal = document.getElementById('notes-modal');
const notesDragHandle = document.getElementById('notes-drag-handle');
const notesMinimizeBtn = document.getElementById('notes-minimize-btn');
const notesCloseBtn = document.getElementById('notes-close-btn');
const notesBodyContent = document.getElementById('notes-body-content');
const projectNotesInput = document.getElementById('project-notes-input');

// Rampas de Temperatura
const tempRampModal = document.getElementById('temp-ramp-modal');
const tempRampDragHandle = document.getElementById('temp-ramp-drag-handle');
const tempRampMinimizeBtn = document.getElementById('temp-ramp-minimize-btn');
const tempRampCloseBtn = document.getElementById('temp-ramp-close-btn');
const tempRampBodyContent = document.getElementById('temp-ramp-body-content');
const rampSelect = document.getElementById('ramp-select');
const btnAddRamp = document.getElementById('btn-add-ramp');
const btnDelRamp = document.getElementById('btn-del-ramp');
const rampNameInput = document.getElementById('ramp-name-input');
const rampVisibleCheckbox = document.getElementById('ramp-visible-checkbox');
const rampStagesBody = document.getElementById('ramp-stages-body');
const btnAddStage = document.getElementById('btn-add-stage');
const rampSummaryInfo = document.getElementById('ramp-summary-info');

// Dupla Linha Tracejada + Tooltip
const crosshairV = document.getElementById('crosshair-v');
const crosshairH = document.getElementById('crosshair-h');
const crosshairTooltip = document.getElementById('crosshair-tooltip');

// Menu de Contexto (Botão Direito)
const contextMenu = document.getElementById('context-menu');
const yModeRaw = document.getElementById('y-mode-raw');
const yModePercent = document.getElementById('y-mode-percent');
const ctxResetZoom = document.getElementById('ctx-reset-zoom');
const ctxResetXY = document.getElementById('ctx-reset-xy');
const ctxLimXMin = document.getElementById('ctx-lim-x-min');
const ctxLimXMax = document.getElementById('ctx-lim-x-max');
const ctxLimYMin = document.getElementById('ctx-lim-y-min');
const ctxLimYMax = document.getElementById('ctx-lim-y-max');
const ctxOffsetY = document.getElementById('ctx-offset-y');

// Modal Arrastável de Integração
const modal = document.getElementById('integration-modal');
const modalDragHandle = document.getElementById('modal-drag-handle');
const modalMinimizeBtn = document.getElementById('modal-minimize-btn');
const modalCloseBtn = document.getElementById('modal-close-btn');
const modalBodyContent = document.getElementById('modal-body-content');
const btnCancelIntegration = document.getElementById('btn-cancel-integration');
const btnApplyIntegration = document.getElementById('btn-apply-integration');

const startSlider = document.getElementById('start-slider');
const endSlider = document.getElementById('end-slider');
const startIdxFeedback = document.getElementById('start-idx-feedback');
const startTimeFeedback = document.getElementById('start-time-feedback');
const endIdxFeedback = document.getElementById('end-idx-feedback');
const endTimeFeedback = document.getElementById('end-time-feedback');

const peakNameInput = document.getElementById('peak-name-input');
const baselineTypeSelect = document.getElementById('baseline-type-select');
const droplineControls = document.getElementById('dropline-controls');
const tangentControls = document.getElementById('tangent-controls');

const glStartSlider = document.getElementById('gl-start-slider');
const glEndSlider = document.getElementById('gl-end-slider');
const anchorSlider = document.getElementById('anchor-slider');

const targetVisible = document.getElementById('target-visible');
const targetCustom = document.getElementById('target-custom');
const customSelectionBox = document.getElementById('custom-selection-box');
const customChecklist = document.getElementById('custom-sample-checklist');
const btnSelectAll = document.getElementById('btn-select-all');
const btnDeselectAll = document.getElementById('btn-deselect-all');

// Gaveta Inferior de Resultados
const bottomPanel = document.getElementById('bottom-panel');
const panelResizer = document.getElementById('panel-resizer');
const panelTabList = document.getElementById('panel-tab-list');
const peaksTableBody = document.getElementById('peaks-table-body');
const thAreaPct = document.getElementById('th-area-pct');
const thNorm = document.getElementById('th-norm');
const panelMinimizeBtn = document.getElementById('panel-minimize-btn');
const panelCloseBtn = document.getElementById('panel-close-btn');
const btnToggleShading = document.getElementById('btn-toggle-shading');
const btnToggleLabels = document.getElementById('btn-toggle-labels');
const btnExportTableCsv = document.getElementById('btn-export-table-csv');

// ── CONFIGURAÇÕES E ESTADOS GLOBAIS ──
const COLOR_PALETTE = [
    '#2563eb', '#dc2626', '#16a34a', '#f59e0b', '#8b5cf6',
    '#06b6d4', '#ea580c', '#ec4899', '#10b981', '#6366f1',
    '#84cc16', '#14b8a6', '#f43f5e', '#a855f7', '#eab308'
];

let samples = [];
let fidStash = null;   // amostras de GC-FID guardadas enquanto o GC-MS usa `samples` para os "sinais" (curvas do EIC / TIC)
let appliedPeaks = [];

// Pico em edição (janela de integração aberta pela tabela). null = integração nova.
let editingPeak = null;

// O fluxo "clicar início e fim no gráfico" serve a duas ações: integrar um pico ou remover picos de uma região.
let pickPurpose = 'integrate';   // o que o fluxo de seleção início/fim vai fazer
const PICK_TEXT = {              // textos do banner por finalidade: [início, fim]
    'integrate': ['Início do pico', 'Fim do pico'],
    'remove': ['Início da região', 'Fim da região'],
    'ms-avg': ['Início do intervalo', 'Fim do intervalo'],   // GC-MS: espectro médio
    'ms-bg': ['Início do fundo', 'Fim do fundo'],            // GC-MS: região de fundo
    'align-anchor': ['Clique no pico âncora', 'Clique no pico âncora']   // alinhamento: um clique
};
let regionRemoval = null;        // { t0, t1, doomed:Set } enquanto a janela "Remover picos da região" está aberta
const selectedPeaks = new Set(); // linhas marcadas na tabela de picos

// Histórico de desfazer/refazer dos picos. Os objetos de pico são tratados como imutáveis
// (qualquer alteração cria um novo objeto), então cada snapshot é só uma cópia rasa do array.
const HISTORY_LIMIT = 100;
let undoStack = [];
let redoStack = [];
let currentTheme = 'light';
let activePanelTab = 'all';

let isYPercentMode = false;
let currentOffsetY = 0;

let isShadingVisible = true;
let isLabelsVisible = true;

let pickState = 'idle';
let pickedStartTime = null;

// Bloco de Notas (atrelado ao projeto)
let projectNotes = '';

// Rampas de Temperatura do Forno GC
let temperatureRamps = [
    {
        id: 'ramp_1',
        name: 'Rampa 1',
        visible: false,
        stages: [
            { rate: null, temp: 40, hold: 2 },
            { rate: 10, temp: 200, hold: 5 }
        ]
    }
];
let selectedRampId = 'ramp_1';

// Nome e caminho do projeto atualmente aberto no diretório
let currentOpenProjectFilename = null;
let currentOpenProjectPath = null;

// ── Modo GC-FID / GC-MS ──
let appMode = 'gcfid';          // 'gcfid' | 'gcms'
let msSamples = [];             // arquivos GC-MS do projeto (apenas um é exibido por vez)
let activeMsId = null;          // id do arquivo GC-MS em exibição
let msYRelative = true;         // espectro em abundância relativa (%) ou absoluta
let msCurrentSpectrum = null;   // último espectro de scan carregado/exibido
let msPanelMode = 'scan';       // 'scan' (um scan) | 'avg' (espectro médio de um intervalo)
let msAvgResult = null;         // resposta de /api/gcms/average para a amostra ativa
let msReference = null;         // espectro de referência { name, mz, intensity, source, meta } para a aba Comparar
let msCmpSource = 'scan';       // o que a aba Comparar compara com a referência: 'scan' | 'avg'
let msCmpResult = null;         // última resposta de /api/gcms/compare
// EIC é global (não por amostra): o mesmo íon é comparado ao trocar de corrida
let msEic = { ions: [], tol: 0.5, showTic: false };
let msEicData = null;           // { key, series } da amostra ativa (séries por scan, uma por íon)
const EIC_MAX = 8;
// Normalização das corridas (global): fator pelo padrão interno ou manual; coluna da tabela e/ou escala de exibição
let normSettings = { mode: 'none', isName: '', metric: 'area', plot: false };
// Método de compostos (lista reutilizável): { id, name, rt, win, baseline, ion, tol, is, domain: 'fid'|'ms' }
let compounds = [];
let cmpSeq = 1;
// Grupos de amostras (organização): { id, name, collapsed, domain: 'fid'|'ms' }; cada amostra guarda `group` (id)
let sampleGroups = [];
const EIC_COLORS = ['#e11d48', '#0284c7', '#16a34a', '#d97706', '#7c3aed', '#0891b2', '#db2777', '#65a30d'];

// Valida o intervalo do espectro médio vindo de um arquivo de projeto
function sanitizeMsAvg(m) {
    const a = m.avg, n = (m.time || []).length;
    if (!a || !(a.i0 >= 0 && a.i1 < n && a.i0 <= a.i1)) return null;
    const bgOk = a.bg && a.bg.i0 >= 0 && a.bg.i1 < n && a.bg.i0 <= a.bg.i1;
    return {
        i0: a.i0, i1: a.i1,
        bg: bgOk ? { i0: a.bg.i0, i1: a.bg.i1 } : null,
        bin: [1, 0.5, 0.1, 0.01].includes(a.bin) ? a.bin : 1,
        subtract: a.subtract !== false
    };
}

// ── Sinais integráveis ──
// O código de integração, tabela de picos e edição trabalha sobre `samples`. No GC-FID são as amostras.
// No GC-MS, `samples` passa a ser a lista de "sinais" da amostra ativa: as curvas do EIC (uma por íon)
// ou, sem EIC, o TIC. Cada sinal tem id "msig|<arquivo>|<íon ou TIC>|<tolerância>", que o servidor entende.
function msSignalId(msId, ion, tol) {
    return `msig|${msId}|${ion == null ? 'TIC' : ion}|${ion == null ? 0 : tol}`;
}

function msSignalName(s, ion) {
    return ion == null ? `${s.name} · TIC` : `${s.name} · m/z ${+ion.toFixed(2)}`;
}

function makeMsSignal(s, ion, tol, color, intensity) {
    return {
        id: msSignalId(s.id, ion, tol), name: msSignalName(s, ion), color,
        time: s.time, intensity, visible: true, maxVal: 0, ms_id: s.id, ion, tol
    };
}

function buildMsSignals() {
    const s = getActiveMs();
    if (!s) return [];
    if (isEicActive()) {
        msEicData.signals.forEach((sig, k) => { sig.name = msSignalName(s, msEic.ions[k]); });
        return msEicData.signals;
    }
    if (!s._ticSig || s._ticSig.time !== s.time) s._ticSig = makeMsSignal(s, null, 0, s.color, s.tic);
    s._ticSig.name = msSignalName(s, null);
    s._ticSig.color = s.color;
    return [s._ticSig];
}

function rebuildMsSignals() {
    if (appMode === 'gcms') samples = buildMsSignals();
}

function useFidDomain() {
    if (fidStash !== null) { samples = fidStash; fidStash = null; }
}

// Alinha `samples` ao modo atual (chamado ao trocar de modo e depois de limpar/abrir projetos)
function activateDomain() {
    if (appMode === 'gcms') {
        if (fidStash === null) fidStash = samples;
        samples = buildMsSignals();
    } else {
        useFidDomain();
    }
}

function getFidSamples() {
    return fidStash !== null ? fidStash : samples;
}

// Picos cujo sinal está na tela agora (os de outras curvas/modos ficam guardados, só não aparecem)
function domainPeaks() {
    return appliedPeaks.filter(pk => samples.some(s => s.id === pk.sample_id));
}

function eicKey(s) {
    return `${s.id}|${msEic.tol}|${msEic.ions.join(',')}`;
}

// EIC ligado e com dados prontos para a amostra ativa (enquanto carrega, o gráfico mostra o TIC)
function isEicActive() {
    const s = getActiveMs();
    return appMode === 'gcms' && !!s && msEic.ions.length > 0 && !!msEicData && msEicData.key === eicKey(s);
}

function getYTitle() {
    if (appMode === 'gcms') {
        const eic = isEicActive();
        return isYPercentMode ? (eic ? 'EIC Relativo (%)' : 'TIC Relativo (%)') : (eic ? 'EIC (abundância)' : 'TIC (abundância)');
    }
    return isYPercentMode ? 'Intensidade Relativa (%)' : 'Intensidade / Sinal FID';
}

function getPlotLayout() {
    const isDark = currentTheme === 'dark';
    const hasVisibleRamp = (temperatureRamps || []).some(r => r.visible);

    const layout = {
        uirevision: 'constant_zoom_lock',
        paper_bgcolor: isDark ? '#090d16' : '#ffffff',
        plot_bgcolor: isDark ? '#090d16' : '#ffffff',
        font: { color: isDark ? '#cbd5e1' : '#334155', size: 12 },
        margin: { l: 65, r: hasVisibleRamp ? 55 : 25, t: 25, b: 50 },
        showlegend: false,
        hovermode: false,
        xaxis: {
            title: 'Tempo de Retenção (min)',
            gridcolor: isDark ? '#1a243b' : '#f1f5f9',
            zerolinecolor: isDark ? '#233152' : '#e2e8f0',
            tickformat: '.2f'
        },
        yaxis: {
            title: getYTitle(),
            gridcolor: isDark ? '#1a243b' : '#f1f5f9',
            zerolinecolor: isDark ? '#233152' : '#e2e8f0'
        }
    };

    if (hasVisibleRamp) {
        layout.yaxis2 = {
            title: { text: 'Temperatura (°C)', font: { color: isDark ? '#94a3b8' : '#64748b', size: 11 } },
            tickfont: { color: isDark ? '#94a3b8' : '#64748b', size: 10 },
            overlaying: 'y',
            side: 'right',
            showgrid: false,
            zeroline: false
        };
    }

    // GC-MS: linha tracejada na posição do scan cujo espectro está na gaveta inferior
    if (isEicActive()) {
        layout.showlegend = true;
        layout.legend = { orientation: 'h', x: 0, y: 1.0, yanchor: 'bottom', font: { size: 11 }, bgcolor: 'rgba(0,0,0,0)' };
        layout.margin.t = 40;
    }

    if (appMode === 'gcms') layout.shapes = getMsMarkerShapes();
    else if (regionRemoval) layout.shapes = getRegionShapes();

    return layout;
}

const plotConfig = {
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['select2d', 'lasso2d']
};

// No GC-MS o duplo clique escolhe o scan (em vez de resetar o zoom; o reset fica no menu de contexto)
function getPlotConfig() {
    return appMode === 'gcms' ? { ...plotConfig, doubleClick: false } : plotConfig;
}

Plotly.newPlot(plotDiv, [], getPlotLayout(), getPlotConfig());

// ── CORREÇÃO DE ESCALA Y ──
function scaleY(val, sample, sampleIndex = 0) {
    if (appMode === 'gcms') sampleIndex = 0;   // sem deslocamento vertical entre as curvas do GC-MS
    let y = val;
    const maxVal = sample.maxVal || arrayMax(sample.intensity) || 1;

    if (isYPercentMode) {
        y = (val / maxVal) * 100.0;
        if (currentOffsetY > 0) {
            const offsetPct = Math.min(50, (currentOffsetY / maxVal) * 100.0);
            y += (sampleIndex * offsetPct);
        }
    } else {
        y = val * (sample.dispFactor || 1);   // normalização das corridas (só exibição)
        if (currentOffsetY > 0) {
            y += (sampleIndex * currentOffsetY);
        }
    }
    return y;
}

// Máximo de um array grande. Math.max(...arr) lança RangeError acima de ~125 mil elementos.
function arrayMax(arr) {
    let m = -Infinity;
    for (let i = 0; i < arr.length; i++) if (arr[i] > m) m = arr[i];
    return m;
}

function formatIntensity(val) {
    if (isYPercentMode) return val.toFixed(2) + ' %';
    const abs = Math.abs(val);
    if (abs >= 1e6) return (val / 1e6).toFixed(2) + 'M';
    if (abs >= 1e3) return (val / 1e3).toFixed(1) + 'k';
    return val.toFixed(1);
}

// ── CÁLCULO E INTERPOLAÇÃO DE RAMPAS DE TEMPERATURA (GC) ──
function calculateRampPoints(ramp, maxPlotTime = 0) {
    const stages = ramp.stages || [];
    if (stages.length === 0) return { x: [], y: [] };

    const x = [];
    const y = [];

    let t = 0;
    let currentTemp = Number(stages[0].temp) || 40;

    // Estágio inicial (Hold inicial a currentTemp)
    x.push(t);
    y.push(currentTemp);

    const initialHold = Number(stages[0].hold) || 0;
    if (initialHold > 0) {
        t += initialHold;
        x.push(t);
        y.push(currentTemp);
    }

    // Estágios subsequentes (aquecimento/resfriamento + hold)
    for (let i = 1; i < stages.length; i++) {
        const st = stages[i];
        const targetTemp = Number(st.temp) || currentTemp;
        const rate = Number(st.rate) || 10;
        const hold = Number(st.hold) || 0;

        if (rate > 0 && Math.abs(targetTemp - currentTemp) > 0.001) {
            const rampDuration = Math.abs(targetTemp - currentTemp) / rate;
            t += rampDuration;
            currentTemp = targetTemp;
            x.push(t);
            y.push(currentTemp);
        } else {
            currentTemp = targetTemp;
        }

        if (hold > 0) {
            t += hold;
            x.push(t);
            y.push(currentTemp);
        }
    }

    const finalRampTime = t;
    const finalRampTemp = currentTemp;

    // Se o cromatograma durar mais que a rampa, estende em linha horizontal
    if (maxPlotTime > t) {
        x.push(maxPlotTime);
        y.push(currentTemp);
    }

    return { x, y, finalTime: finalRampTime, finalTemp: finalRampTemp };
}

function calculateTemperatureAtTime(ramp, tMin) {
    const pts = calculateRampPoints(ramp, 0);
    const x = pts.x;
    const y = pts.y;
    if (!x || x.length === 0) return null;

    if (tMin <= x[0]) return y[0];
    if (tMin >= x[x.length - 1]) return y[y.length - 1];

    for (let i = 0; i < x.length - 1; i++) {
        if (tMin >= x[i] && tMin <= x[i + 1]) {
            const spanX = x[i + 1] - x[i];
            if (spanX <= 0.00001) return y[i + 1];
            const frac = (tMin - x[i]) / spanX;
            return y[i] + (y[i + 1] - y[i]) * frac;
        }
    }
    return y[y.length - 1];
}

// ── RETICULADO: TEMPO APENAS EM MINUTOS + INTENSIDADE + TEMPERATURA ──
chartArea.addEventListener('mousemove', (e) => {
    const fullLayout = plotDiv._fullLayout;
    if (!fullLayout || !fullLayout.xaxis || !fullLayout.yaxis) return;

    const rect = plotDiv.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;
    const xaxis = fullLayout.xaxis;
    const yaxis = fullLayout.yaxis;

    const inPlotX = mouseX >= xaxis._offset && mouseX <= xaxis._offset + xaxis._length;
    const inPlotY = mouseY >= yaxis._offset && mouseY <= yaxis._offset + yaxis._length;

    if (!inPlotX || !inPlotY) {
        crosshairV.style.display = 'none';
        crosshairH.style.display = 'none';
        crosshairTooltip.style.display = 'none';
        return;
    }

    crosshairV.style.display = 'block';
    crosshairH.style.display = 'block';
    crosshairV.style.left = `${mouseX}px`;
    crosshairH.style.top = `${mouseY}px`;

    const timeMin = xaxis.p2c(mouseX - xaxis._offset);
    const intensity = yaxis.p2c(mouseY - yaxis._offset);

    let tempTooltipHtml = '';
    const visibleRamps = (temperatureRamps || []).filter(r => r.visible);
    if (visibleRamps.length > 0) {
        visibleRamps.forEach(r => {
            const tempVal = calculateTemperatureAtTime(r, timeMin);
            if (tempVal !== null) {
                const label = visibleRamps.length > 1 ? `Temp (${r.name}):` : 'Temperatura:';
                tempTooltipHtml += `<br>${label} ${tempVal.toFixed(1)} °C`;
            }
        });
    }

    crosshairTooltip.style.display = 'block';
    crosshairTooltip.innerHTML = `
        Tempo: ${timeMin.toFixed(3)} min<br>
        Intensidade: ${formatIntensity(intensity)}${tempTooltipHtml}
    `;

    let tooltipLeft = mouseX + 16;
    let tooltipTop = mouseY + 16;
    if (tooltipLeft + 200 > rect.width) tooltipLeft = mouseX - 210;
    if (tooltipTop + 75 > rect.height) tooltipTop = mouseY - 80;

    crosshairTooltip.style.left = `${tooltipLeft}px`;
    crosshairTooltip.style.top = `${tooltipTop}px`;
});

chartArea.addEventListener('mouseleave', () => {
    crosshairV.style.display = 'none';
    crosshairH.style.display = 'none';
    crosshairTooltip.style.display = 'none';
});

// Scroll de Intensidade Y
plotDiv.addEventListener('wheel', (e) => {
    e.preventDefault();
    const fullLayout = plotDiv._fullLayout;
    if (!fullLayout || !fullLayout.yaxis) return;

    const yaxis = fullLayout.yaxis;
    const currentY = yaxis.range.slice();
    let yMin = Math.max(0, currentY[0]);
    let yMax = currentY[1];
    const span = yMax - yMin;

    const zoomFactor = e.deltaY < 0 ? 0.82 : 1.22;
    const newSpan = Math.max(isYPercentMode ? 1 : 10, span * zoomFactor);

    Plotly.relayout(plotDiv, {
        'yaxis.range': [0, newSpan],
        'yaxis.autorange': false
    });
}, { passive: false });

function sanitizeFilename(filename) {
    return filename.replace(/\.cdf$/i, '');
}

function getTargetSamples() {
    if (editingPeak) return samples.filter(s => s.id === editingPeak.sample_id);
    if (targetVisible.checked) {
        return samples.filter(s => s.visible);
    } else {
        const checkedBoxes = customChecklist.querySelectorAll('input:checked');
        const checkedIds = Array.from(checkedBoxes).map(cb => cb.value);
        return samples.filter(s => checkedIds.includes(s.id));
    }
}

// Traços dos picos já integrados (polígono, linha de base, ápice e rótulo) dos sinais visíveis.
// Compartilhado entre GC-FID e GC-MS (EIC/TIC).
function getAppliedPeakTraces() {
    const plotData = [];
    const isDark = currentTheme === 'dark';
    appliedPeaks.forEach((pk) => {
        if (pk === editingPeak) return;   // o preview da edição já o redesenha
        const parent = samples.find(s => s.id === pk.sample_id);
        if (!parent || !parent.visible) return;

        const sIdx = samples.indexOf(parent);
        const runPeaks = appliedPeaks.filter(p => p.sample_id === pk.sample_id);
        const peakNum = runPeaks.indexOf(pk) + 1;
        const peakLabel = `#${peakNum} ${pk.peak_name}`;

        const scaledPolyY = pk.polygon.y.map(v => scaleY(v, parent, sIdx));
        const scaledBaseY = pk.baseline_line.y.map(v => scaleY(v, parent, sIdx));
        const scaledApexY = scaleY(pk.apex.y, parent, sIdx);
        const doomed = !!(regionRemoval && regionRemoval.doomed.has(pk));   // será removido: destaque em vermelho

        if (isShadingVisible) {
            plotData.push({
                x: pk.polygon.x,
                y: scaledPolyY,
                mode: 'lines',
                fill: 'toself',
                fillcolor: doomed ? '#ef444477' : parent.color + '40',
                line: { width: 0 },
                hoverinfo: 'none'
            });

            plotData.push({
                x: pk.baseline_line.x,
                y: scaledBaseY,
                mode: 'lines',
                line: { color: parent.color, dash: 'dash', width: 1.5 },
                hoverinfo: 'none'
            });
        }

        plotData.push({
            x: [pk.apex.x],
            y: [scaledApexY],
            mode: isLabelsVisible ? 'markers+text' : 'markers',
            marker: { size: doomed ? 9 : 7, color: doomed ? '#ef4444' : parent.color },
            text: isLabelsVisible ? [peakLabel] : undefined,
            textposition: 'top center',
            textfont: {
                size: 11,
                color: isDark ? '#f8fafc' : '#0f172a',
                family: 'sans-serif'
            },
            hoverinfo: 'none'
        });
    });
    return plotData;
}

// Prévia da integração (janela aberta) e marcação do ponto inicial durante a seleção.
function getIntegrationPreviewTraces() {
    const plotData = [];
    if (modal.style.display !== 'none' && samples.length > 0) {
        updateOverlapWarning();
        const targetList = getTargetSamples();
        const baseType = baselineTypeSelect.value;

        const sIdxRaw = parseInt(startSlider.value);
        const eIdxRaw = parseInt(endSlider.value);
        const glSRaw = parseInt(glStartSlider.value);
        const glERaw = parseInt(glEndSlider.value);
        const ancRaw = parseInt(anchorSlider.value);

        targetList.forEach(s => {
            const sIdxSample = samples.indexOf(s);
            const n = s.time.length;
            const sIdx = Math.max(0, Math.min(sIdxRaw, n - 1));
            const eIdx = Math.max(0, Math.min(eIdxRaw, n - 1));
            if (sIdx >= eIdx) return;

            const xSeg = s.time.slice(sIdx, eIdx + 1);
            const ySegRaw = s.intensity.slice(sIdx, eIdx + 1);
            const ySeg = ySegRaw.map(v => scaleY(v, s, sIdxSample));

            let yBase = [];
            let blLineX = [];
            let blLineY = [];

            if (baseType === 'Valley-to-Valley') {
                const y0 = ySeg[0];
                const y1 = ySeg[ySeg.length - 1];
                for (let i = 0; i < xSeg.length; i++) {
                    yBase.push(y0 + (y1 - y0) * (i / (xSeg.length - 1)));
                }
                blLineX = [xSeg[0], xSeg[xSeg.length - 1]];
                blLineY = [y0, y1];

                plotData.push({
                    x: [xSeg[0], xSeg[xSeg.length - 1]],
                    y: [y0, y1],
                    mode: 'markers',
                    marker: { symbol: 'x', size: 10, color: s.color, line: { width: 2 } },
                    hoverinfo: 'none'
                });

            } else if (baseType === 'Drop Line') {
                const glS = Math.max(0, Math.min(glSRaw, n - 1));
                const glE = Math.max(0, Math.min(glERaw, n - 1));
                const glX0 = s.time[glS], glY0 = scaleY(s.intensity[glS], s, sIdxSample);
                const glX1 = s.time[glE], glY1 = scaleY(s.intensity[glE], s, sIdxSample);

                blLineX = [glX0, glX1];
                blLineY = [glY0, glY1];

                for (let i = 0; i < xSeg.length; i++) {
                    const frac = (xSeg[i] - glX0) / (glX1 - glX0 || 1);
                    yBase.push(glY0 + (glY1 - glY0) * frac);
                }

                plotData.push({
                    x: [glX0, glX1],
                    y: [glY0, glY1],
                    mode: 'markers',
                    marker: { symbol: 'circle', size: 8, color: s.color },
                    hoverinfo: 'none'
                });

                plotData.push({
                    x: [xSeg[0], xSeg[xSeg.length - 1]],
                    y: [ySeg[0], ySeg[ySeg.length - 1]],
                    mode: 'markers',
                    marker: { symbol: 'x', size: 10, color: s.color, line: { width: 2 } },
                    hoverinfo: 'none'
                });

                plotData.push({
                    x: [xSeg[0], xSeg[0]],
                    y: [ySeg[0], yBase[0]],
                    mode: 'lines',
                    line: { color: s.color, dash: 'dot', width: 1.5 },
                    hoverinfo: 'none'
                });
                plotData.push({
                    x: [xSeg[xSeg.length - 1], xSeg[xSeg.length - 1]],
                    y: [ySeg[ySeg.length - 1], yBase[yBase.length - 1]],
                    mode: 'lines',
                    line: { color: s.color, dash: 'dot', width: 1.5 },
                    hoverinfo: 'none'
                });

            } else if (baseType === 'Tangent Skim') {
                const anc = Math.max(0, Math.min(ancRaw, n - 1));
                const ancX = s.time[anc], ancY = scaleY(s.intensity[anc], s, sIdxSample);
                const startX = xSeg[0], startY = ySeg[0];

                blLineX = [startX, ancX];
                blLineY = [startY, ancY];

                for (let i = 0; i < xSeg.length; i++) {
                    const frac = (xSeg[i] - startX) / (ancX - startX || 1);
                    yBase.push(startY + (ancY - startY) * frac);
                }

                plotData.push({
                    x: [ancX],
                    y: [ancY],
                    mode: 'markers',
                    marker: { symbol: 'circle', size: 8, color: s.color },
                    hoverinfo: 'none'
                });

                plotData.push({
                    x: [startX],
                    y: [startY],
                    mode: 'markers',
                    marker: { symbol: 'x', size: 10, color: s.color, line: { width: 2 } },
                    hoverinfo: 'none'
                });
            }

            plotData.push({
                x: blLineX,
                y: blLineY,
                mode: 'lines',
                line: { color: s.color, dash: 'dash', width: 1.5 },
                hoverinfo: 'none'
            });

            const polyX = xSeg.concat(xSeg.slice().reverse());
            const polyY = ySeg.concat(yBase.slice().reverse());

            plotData.push({
                x: polyX,
                y: polyY,
                mode: 'lines',
                fill: 'toself',
                fillcolor: s.color + '45',
                line: { width: 0 },
                hoverinfo: 'none'
            });
        });
    }

    if (pickState === 'end' && pickedStartTime !== null && (pickPurpose === 'integrate' || pickPurpose === 'remove')) {
        samples.filter(s => s.visible).forEach(s => {
            const sIdxSample = samples.indexOf(s);
            const idx = findClosestIndex(s, pickedStartTime);
            plotData.push({
                x: [s.time[idx]],
                y: [scaleY(s.intensity[idx], s, sIdxSample)],
                mode: 'markers',
                marker: { symbol: 'x', size: 11, color: s.color, line: { width: 2.5 } },
                hoverinfo: 'none'
            });
        });
    }
    return plotData;
}

// ── Atualização do Gráfico ──
function updateApp() {
    if (appMode === 'gcms') { updateMsApp(); return; }
    sampleCount.innerText = samples.length;
    const plotData = [];
    const isDark = currentTheme === 'dark';

    applyDisplayFactors();
    samples.forEach((s, sIdx) => {
        s.maxVal = arrayMax(s.intensity);
        const scaledY = s.intensity.map(v => scaleY(v, s, sIdx));

        plotData.push({
            x: s.time,
            y: scaledY,
            name: s.name,
            type: 'scattergl',
            mode: 'lines',
            line: { color: s.color, width: 1.2 },
            visible: s.visible
        });
    });

    getAppliedPeakTraces().forEach(tr => plotData.push(tr));
    getIntegrationPreviewTraces().forEach(tr => plotData.push(tr));

    // Adiciona as rampas de temperatura visíveis
    const visibleRamps = (temperatureRamps || []).filter(r => r.visible);
    if (visibleRamps.length > 0) {
        let maxPlotTime = 0;
        samples.filter(s => s.visible).forEach(s => {
            if (s.time && s.time.length > 0) {
                const lastT = s.time[s.time.length - 1];
                if (lastT > maxPlotTime) maxPlotTime = lastT;
            }
        });

        visibleRamps.forEach(ramp => {
            const pts = calculateRampPoints(ramp, maxPlotTime);
            if (pts.x && pts.x.length > 0) {
                plotData.push({
                    x: pts.x,
                    y: pts.y,
                    name: ramp.name || 'Rampa de Temp',
                    type: 'scatter',
                    mode: 'lines',
                    line: {
                        dash: 'dot',
                        width: 1.8,
                        color: isDark ? '#64748b' : '#94a3b8'
                    },
                    yaxis: 'y2',
                    hoverinfo: 'none'
                });
            }
        });
    }

    Plotly.react(plotDiv, plotData, getPlotLayout(), getPlotConfig());
    renderSampleList();
}

function findClosestIndex(sample, targetTimeMin) {
    let closest = 0, minDiff = Infinity;
    for (let i = 0; i < sample.time.length; i++) {
        const diff = Math.abs(sample.time[i] - targetTimeMin);
        if (diff < minDiff) {
            minDiff = diff;
            closest = i;
        }
    }
    return closest;
}

function renderSampleList() {
    if (appMode === 'gcms') { renderMsSampleList(); return; }
    renderGroupedList(samples, 'fid', (s, inGroup) => {
        const li = document.createElement('li');
        li.className = `sample-item ${s.visible ? '' : 'is-hidden'} ${inGroup ? 'in-group' : ''}`;

        li.innerHTML = `
            <div class="sample-left">
                <div class="color-picker-wrapper" title="Mudar cor">
                    <span class="color-dot" style="background-color: ${s.color}"></span>
                    <input type="color" class="color-input-hidden" value="${s.color}">
                </div>
                <span class="sample-name" title="Duplo clique p/ renomear">${escapeHtml(s.name)}</span>
            </div>
            <button class="btn-close" title="Remover amostra">✕</button>
        `;

        // Ação do clique adiada: se virar duplo clique (renomear), é cancelada
        // (senão o 1º clique redesenha a lista e o duplo clique nunca chega ao nome).
        let clickTimer = null;
        li.addEventListener('click', (e) => {
            if (e.target.closest('.color-picker-wrapper') || e.target.closest('.rename-input')) return;
            if (e.detail > 1) return;
            clearTimeout(clickTimer);
            clickTimer = setTimeout(() => {
                s.visible = !s.visible;
                updateApp();
                renderBottomPanel();
            }, 230);
        });

        const colorWrapper = li.querySelector('.color-picker-wrapper');
        const colorInput = li.querySelector('.color-input-hidden');
        colorWrapper.addEventListener('click', (e) => { e.stopPropagation(); colorInput.click(); });
        colorInput.addEventListener('change', (e) => { s.color = e.target.value; updateApp(); renderBottomPanel(); });

        const nameSpan = li.querySelector('.sample-name');
        nameSpan.addEventListener('dblclick', (e) => {
            e.stopPropagation();
            clearTimeout(clickTimer);
            const currentName = s.name;
            const input = document.createElement('input');
            input.type = 'text';
            input.className = 'rename-input';
            input.value = currentName;
            nameSpan.replaceWith(input);
            input.focus();
            input.select();
            const finish = () => { s.name = input.value.trim() || currentName; updateApp(); renderBottomPanel(); };
            input.addEventListener('keydown', (ke) => { if (ke.key === 'Enter') finish(); if (ke.key === 'Escape') updateApp(); });
            input.addEventListener('blur', finish);
            li.draggable = false;   // evita arrastar ao selecionar o texto
        });

        const closeBtn = li.querySelector('.btn-close');
        closeBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            samples = samples.filter(item => item.id !== s.id);
            appliedPeaks = appliedPeaks.filter(item => item.sample_id !== s.id);
            resetPeakState();   // o conjunto de amostras mudou: snapshots antigos ficariam inconsistentes
            updateApp();
            renderBottomPanel();
        });

        return li;
    });
}

// ── MODO INTEGRAÇÃO & ATALHO DE TECLADO 'F' ──
function toggleIntegrationMode() {
    if (samples.length === 0) {
        alert('Carregue arquivos .CDF antes de integrar.');
        return;
    }
    togglePicking('integrate');
}

function toggleRegionRemoveMode() {
    if (domainPeaks().length === 0) {
        showToast('Não há picos integrados para remover.');
        return;
    }
    togglePicking('remove');
}

// Clicar de novo no mesmo botão cancela; clicar no outro troca a finalidade
function togglePicking(purpose) {
    if (pickState !== 'idle' && pickPurpose === purpose) { cancelPicking(); return; }
    if (pickState !== 'idle') cancelPicking();
    startPicking(purpose);
}

integrateBtn.addEventListener('click', toggleIntegrationMode);

// ATALHOS GLOBAIS (F para integrar, Ctrl+S para salvar)
document.addEventListener('keydown', (e) => {
    // Atalho Ctrl+S (Salvar) e Ctrl+Shift+S (Salvar Como)
    if ((e.ctrlKey || e.metaKey) && !e.altKey && (e.code === 'KeyS' || e.key.toLowerCase() === 's')) {
        e.preventDefault();
        if (e.repeat) return;   // segurar a tecla não deve abrir vários diálogos
        if (e.shiftKey) {
            btnSaveAsProject.click();
        } else {
            btnSaveProject.click();
        }
        return;
    }

    const tag = (e.target.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'textarea' || e.target.isContentEditable) {
        return;
    }

    if (e.key === 'f' || e.key === 'F') {
        e.preventDefault();
        toggleIntegrationMode();
    }

    if (e.key === 'Escape') {
        if (pickState !== 'idle') cancelPicking();
        if (modal.style.display !== 'none') closeIntegrationModal();
        closeRemoveRegionModal();
        closeEicModal();
        closeAlignModal();
        closeCompoundsModal();
        closeSpecrefModal();
    }
});

function startPicking(purpose = 'integrate') {
    pickPurpose = purpose;
    pickState = 'start';
    pickedStartTime = null;

    ({ 'remove': regionBtn, 'ms-avg': msAvgBtn, 'ms-bg': msBgBtn, 'align-anchor': alignBtn }[purpose] || integrateBtn).classList.add('active-mode');
    plotDiv.classList.add('picking-mode');
    integrationBanner.style.display = 'flex';
    bannerText.innerText = PICK_TEXT[purpose][0];
    modal.style.display = 'none';
    setIntegrationModalMode(false);
    closeRemoveRegionModal();
}

function cancelPicking() {
    pickState = 'idle';
    integrateBtn.classList.remove('active-mode');
    regionBtn.classList.remove('active-mode');
    msAvgBtn.classList.remove('active-mode');
    msBgBtn.classList.remove('active-mode');
    alignBtn.classList.remove('active-mode');
    plotDiv.classList.remove('picking-mode');
    integrationBanner.style.display = 'none';
    updateApp();
}

bannerCancelBtn.addEventListener('click', cancelPicking);

plotDiv.addEventListener('click', (e) => {
    if (pickState === 'idle') return;

    const fullLayout = plotDiv._fullLayout;
    if (!fullLayout || !fullLayout.xaxis) return;

    const rect = plotDiv.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const xaxis = fullLayout.xaxis;

    if (mouseX < xaxis._offset || mouseX > xaxis._offset + xaxis._length) return;

    const clickedTimeMin = xaxis.p2c(mouseX - xaxis._offset);

    if (pickState === 'start' && pickPurpose === 'align-anchor') {
        cancelPicking();
        alignAnchor.value = clickedTimeMin.toFixed(3);
        detectAlignment();
        return;
    }

    if (pickState === 'start') {
        pickedStartTime = clickedTimeMin;
        pickState = 'end';
        bannerText.innerText = PICK_TEXT[pickPurpose][1];
        updateApp();
    } else if (pickState === 'end') {
        const timeStart = Math.min(pickedStartTime, clickedTimeMin);
        const timeEnd = Math.max(pickedStartTime, clickedTimeMin);

        if (pickPurpose === 'remove') {
            cancelPicking();
            openRemoveRegionModal(timeStart, timeEnd);
            return;
        }

        if (pickPurpose === 'ms-avg' || pickPurpose === 'ms-bg') {
            const purpose = pickPurpose;
            cancelPicking();
            applyMsInterval(purpose, timeStart, timeEnd);
            return;
        }

        const ref = samples.find(s => s.visible) || samples[0];
        const sIdx = findClosestIndex(ref, timeStart);
        const eIdx = findClosestIndex(ref, timeEnd);

        cancelPicking();
        openIntegrationModal(sIdx, eIdx);
    }
});

// ── MENU DE CONTEXTO (BOTÃO DIREITO NO GRÁFICO) ──
plotDiv.addEventListener('contextmenu', (e) => {
    e.preventDefault();
    const rect = chartArea.getBoundingClientRect();
    const posX = e.clientX - rect.left;
    const posY = e.clientY - rect.top;

    contextMenu.style.left = `${posX}px`;
    contextMenu.style.top = `${posY}px`;
    contextMenu.style.display = 'flex';
});

document.addEventListener('click', (e) => {
    if (!contextMenu.contains(e.target)) {
        contextMenu.style.display = 'none';
    }
});

yModeRaw.addEventListener('change', () => {
    isYPercentMode = false;
    contextMenu.style.display = 'none';
    Plotly.relayout(plotDiv, {
        'yaxis.autorange': true,
        'yaxis.title': getYTitle()
    });
    updateApp();
});

yModePercent.addEventListener('change', () => {
    isYPercentMode = true;
    contextMenu.style.display = 'none';
    Plotly.relayout(plotDiv, {
        'yaxis.range': [0, 105],
        'yaxis.autorange': false,
        'yaxis.title': getYTitle()
    });
    updateApp();
});

ctxResetZoom.addEventListener('click', () => {
    contextMenu.style.display = 'none';
    Plotly.relayout(plotDiv, { 'xaxis.autorange': true, 'yaxis.autorange': true });
});

ctxResetXY.addEventListener('click', () => {
    ctxLimXMin.value = '';
    ctxLimXMax.value = '';
    ctxLimYMin.value = '';
    ctxLimYMax.value = '';
    Plotly.relayout(plotDiv, { 'xaxis.autorange': true, 'yaxis.autorange': true });
});

function applyManualLimits() {
    const x0 = parseFloat(ctxLimXMin.value);
    const x1 = parseFloat(ctxLimXMax.value);
    const y0 = parseFloat(ctxLimYMin.value);
    const y1 = parseFloat(ctxLimYMax.value);

    const updates = {};
    if (!isNaN(x0) && !isNaN(x1)) updates['xaxis.range'] = [x0, x1];
    if (!isNaN(y0) && !isNaN(y1)) updates['yaxis.range'] = [y0, y1];

    if (Object.keys(updates).length > 0) {
        updates['xaxis.autorange'] = false;
        updates['yaxis.autorange'] = false;
        Plotly.relayout(plotDiv, updates);
    }
}

[ctxLimXMin, ctxLimXMax, ctxLimYMin, ctxLimYMax].forEach(inp => {
    inp.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') applyManualLimits();
    });
});

ctxOffsetY.addEventListener('input', () => {
    currentOffsetY = parseFloat(ctxOffsetY.value) || 0;
    updateApp();
});

// ── DIRETÓRIO DE PROJETOS (PAINEL DIREITO) ──
function updateCurrentProjectDisplay() {
    const titleBase = 'Chromatogram Viewer';
    if (currentOpenProjectFilename) {
        const pName = currentOpenProjectFilename.replace(/\.json$/i, '');
        document.title = `${pName} - ${titleBase}`;
    } else {
        document.title = titleBase;
    }

    if (!currentProjectDisplay) return;
    if (currentOpenProjectFilename) {
        currentProjectDisplay.innerText = currentOpenProjectFilename.replace(/\.json$/i, '');
        currentProjectDisplay.title = currentOpenProjectPath || currentOpenProjectFilename;
        currentProjectDisplay.style.color = '#0284c7';
    } else {
        currentProjectDisplay.innerText = 'Nenhum';
        currentProjectDisplay.title = 'Nenhum projeto ativo';
        currentProjectDisplay.style.color = 'var(--text-muted)';
    }
}

toggleProjectsPanelBtn.addEventListener('click', () => {
    const isCol = rightSidebar.classList.toggle('collapsed');
    rightResizer.classList.toggle('hidden', isCol);
    toggleProjectsPanelBtn.style.display = isCol ? 'flex' : 'none';
    setTimeout(() => Plotly.Plots.resize(plotDiv), 180);
});

closeProjectsBtn.addEventListener('click', () => {
    rightSidebar.classList.add('collapsed');
    rightResizer.classList.add('hidden');
    toggleProjectsPanelBtn.style.display = 'flex';
    setTimeout(() => Plotly.Plots.resize(plotDiv), 180);
});

async function fetchProjectsInfo() {
    try {
        const res = await fetch('/api/projects/info');
        const data = await res.json();
        currentDirDisplay.innerText = data.directory;
        currentDirDisplay.title = data.directory;
        projectsCount.innerText = data.projects.length;

        projectsList.innerHTML = '';
        if (data.projects.length === 0) {
            projectsList.innerHTML = '<li style="font-size:0.75rem; color:var(--text-muted); padding:8px;">Nenhum projeto salvo nesta pasta.</li>';
            return;
        }

        data.projects.forEach(p => {
            const li = document.createElement('li');
            const isActive = currentOpenProjectFilename === p.filename;
            li.className = `project-card ${isActive ? 'active' : ''}`;
            li.innerHTML = `
                <div class="project-info" title="Clique para carregar ${escapeHtml(p.name)}">
                  <span class="project-name">📄 ${escapeHtml(p.name)}</span>
                  <span class="project-meta">${p.modified} • ${p.size_kb} KB</span>
                </div>
                <button class="btn-close" title="Excluir do disco">✕</button>
            `;

            li.querySelector('.project-info').addEventListener('click', () => loadProjectByName(p.filename));
            li.querySelector('.btn-close').addEventListener('click', async (e) => {
                e.stopPropagation();
                if (confirm(`Deseja excluir "${p.name}"?`)) {
                    await fetch(`/api/projects/delete/${encodeURIComponent(p.filename)}`, { method: 'DELETE' });
                    if (currentOpenProjectFilename === p.filename) {
                        currentOpenProjectFilename = null;
                        currentOpenProjectPath = null;
                        updateCurrentProjectDisplay();
                    }
                    fetchProjectsInfo();
                }
            });

            projectsList.appendChild(li);
        });
    } catch (err) {
        console.error('Falha ao obter projetos:', err);
    }
}

// Alterar pasta via janela nativa do Windows Explorer
btnChangeDirectory.addEventListener('click', async () => {
    btnChangeDirectory.innerText = 'Abrindo...';
    try {
        const res = await fetch('/api/projects/pick-directory', { method: 'POST' });
        const data = await res.json();
        if (data.status === 'OK') {
            fetchProjectsInfo();
        } else if (data.status === 'ERROR') {
            throw new Error(data.message || 'Falha ao abrir a janela de seleção.');
        }
    } catch (err) {
        const newDir = prompt('Não foi possível abrir a janela de seleção. Digite o caminho da pasta:', currentDirDisplay.innerText);
        if (newDir && newDir.trim()) {
            await fetch('/api/projects/set-directory', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ directory: newDir.trim() })
            });
            fetchProjectsInfo();
        }
    } finally {
        btnChangeDirectory.innerText = 'Alterar Pasta';
    }
});

// Constrói payload do projeto completo
function buildProjectPayload() {
    return {
        app: "ChromatoView",
        version: "1.0",
        exportDate: new Date().toISOString(),
        currentTheme: currentTheme,
        isYPercentMode: isYPercentMode,
        currentOffsetY: currentOffsetY,
        notes: projectNotes,
        mode: appMode,
        activeMsId: activeMsId,
        msPanelMode: msPanelMode,
        msEic: { ions: msEic.ions, tol: msEic.tol, showTic: msEic.showTic },
        msReference: msReference,
        // GC-MS: salva apenas o TIC; os espectros ficam em cache_gcms/ no servidor
        msSamples: msSamples.map(m => ({
            id: m.id,
            name: m.name,
            color: m.color,
            n_scans: m.n_scans,
            mz_min: m.mz_min,
            mz_max: m.mz_max,
            scanIdx: m.scanIdx,
            shift: m.shift || 0,
            group: m.group || null,
            normManual: m.normManual || 1,
            avg: m.avg || null,   // intervalo do espectro médio, fundo, bin e subtração (o espectro é recalculado ao abrir)
            time: m.time,
            tic: m.tic
        })),
        temperatureRamps: temperatureRamps,
        samples: getFidSamples().map(s => ({
            id: s.id,
            name: s.name,
            color: s.color,
            visible: s.visible,
            n_points: s.n_points || s.time.length,
            sampling_interval_sec: s.sampling_interval_sec || 0.02,
            delay_time_sec: s.delay_time_sec || 0.0,
            shift: s.shift || 0,
            group: s.group || null,
            normManual: s.normManual || 1,
            time: s.time,
            intensity: s.intensity
        })),
        normSettings: { ...normSettings },
        compounds: compounds,
        sampleGroups: sampleGroups,
        appliedPeaks: appliedPeaks
    };
}

// SALVAR: Se tem projeto aberto, salva direto nele. Se NÃO tem, funciona como Salvar Como (perguntando local e nome genérico)
let isSavingNow = false;   // evita salvar/abrir "Salvar como" duas vezes ao mesmo tempo (ex.: Ctrl+S repetido)

btnSaveProject.addEventListener('click', async () => {
    if (isSavingNow) { showToast('Já existe uma janela de salvar aberta (confira a barra de tarefas).'); return; }
    isSavingNow = true;
    try {
        if (!currentOpenProjectFilename && !currentOpenProjectPath) {
            await triggerSaveAs();
            return;
        }
        const projData = buildProjectPayload();
        await saveCurrentProject(projData);
    } finally {
        isSavingNow = false;
    }
});

// SALVAR COMO: Pergunta o local (janela nativa do Windows) e sugere nome genérico
btnSaveAsProject.addEventListener('click', async () => {
    if (isSavingNow) { showToast('Já existe uma janela de salvar aberta (confira a barra de tarefas).'); return; }
    isSavingNow = true;
    try {
        await triggerSaveAs();
    } finally {
        isSavingNow = false;
    }
});

async function triggerSaveAs() {
    const defaultName = currentOpenProjectFilename
        ? currentOpenProjectFilename.replace(/\.json$/i, '') + '_copia'
        : `Projeto_${new Date().toISOString().slice(0, 10)}`;

    try {
        showToast('Abrindo "Salvar como"…');   // feedback imediato (a janela do Windows pode demorar ou ficar atrás do navegador)
        // Abre diálogo nativo do Windows Explorer para escolher pasta e nome
        const res = await fetch('/api/projects/pick-save-path', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ default_name: defaultName })
        });
        const data = await res.json();

        if (data.status === 'OK' && data.filepath) {
            const projData = buildProjectPayload();
            await saveProjectToServer({
                filepath: data.filepath,
                name: data.filename,
                project_data: projData
            });
            return;
        } else if (data.status === 'CANCELLED') {
            showToast('Salvamento cancelado');
            return;
        }
    } catch (err) {
        console.warn('Diálogo nativo de salvar falhou, usando fallback:', err);
    }

    // Fallback: prompt no navegador caso o diálogo nativo não esteja disponível
    const projName = prompt('Salvar Como (digite o nome do projeto):', defaultName);
    if (!projName || !projName.trim()) return;

    const projData = buildProjectPayload();
    await saveProjectToServer({
        name: projName.trim(),
        project_data: projData
    });
}

async function saveCurrentProject(projData) {
    await saveProjectToServer({
        filepath: currentOpenProjectPath,
        name: currentOpenProjectFilename ? currentOpenProjectFilename.replace(/\.json$/i, '') : null,
        project_data: projData
    });
}

async function saveProjectToServer(options) {
    try {
        const res = await fetch('/api/projects/save', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(options)
        });
        const data = await res.json();
        if (data.status === 'OK') {
            currentOpenProjectFilename = data.filename;
            currentOpenProjectPath = data.path;
            updateCurrentProjectDisplay();
            await fetchProjectsInfo();
            const displayName = data.filename.replace(/\.json$/i, '');
            showToast(`Projeto "${displayName}" salvo`);   // não bloqueia como o alert (importante para Ctrl+S)
        } else {
            throw new Error(data.detail || 'Falha na resposta do servidor');
        }
    } catch (err) {
        alert('Erro ao salvar projeto: ' + err.message);
    }
}

// IMPORTAR (Arquivo .json do computador)
btnLoadExternalProject.addEventListener('click', () => {
    externalProjectInput.click();
});

externalProjectInput.addEventListener('change', async (e) => {
    const file = e.target.files && e.target.files[0];
    if (!file) return;

    try {
        const text = await file.text();
        const data = JSON.parse(text);

        if (!data.samples || !Array.isArray(data.samples)) {
            throw new Error('O arquivo selecionado não contém um formato de projeto válido.');
        }

        currentOpenProjectFilename = file.name;
        await applyLoadedProjectData(data);
        updateCurrentProjectDisplay();
        await fetchProjectsInfo();
    } catch (err) {
        alert('Erro ao importar arquivo: ' + err.message);
    } finally {
        externalProjectInput.value = '';
    }
});

// FECHAR PROJETO
btnCloseProject.addEventListener('click', closeCurrentProject);

function closeCurrentProject() {
    if (getFidSamples().length === 0 && appliedPeaks.length === 0 && msSamples.length === 0 && !currentOpenProjectFilename && !currentOpenProjectPath) {
        return;
    }

    const confirmClose = confirm('Deseja fechar o projeto atual? Alterações não salvas serão perdidas.');
    if (!confirmClose) return;

    // Reseta estado local
    useFidDomain();
    samples = [];
    appliedPeaks = [];
    resetPeakState();
    clearMsState();
    activateDomain();
    normSettings = { mode: 'none', isName: '', metric: 'area', plot: false };
    compounds = [];
    sampleGroups = [];
    currentOpenProjectFilename = null;
    currentOpenProjectPath = null;
    currentOffsetY = 0;
    isYPercentMode = false;
    if (yModeRaw) yModeRaw.checked = true;
    if (yModePercent) yModePercent.checked = false;
    if (ctxOffsetY) ctxOffsetY.value = 0;

    // Reseta notas e rampas
    projectNotes = '';
    if (projectNotesInput) projectNotesInput.value = '';
    temperatureRamps = [
        {
            id: 'ramp_1',
            name: 'Rampa 1',
            visible: false,
            stages: [
                { rate: null, temp: 40, hold: 2 },
                { rate: 10, temp: 200, hold: 5 }
            ]
        }
    ];
    selectedRampId = 'ramp_1';
    renderRampModal();
    if (notesModal) notesModal.style.display = 'none';
    if (tempRampModal) tempRampModal.style.display = 'none';
    if (notesBtn) notesBtn.classList.remove('active-mode');
    if (tempRampBtn) tempRampBtn.classList.remove('active-mode');

    // Notifica backend para limpar cache
    fetch('/api/project/sync', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify([])
    }).catch(() => { });

    // Fecha janelas flutuantes
    cancelPicking();
    closeIntegrationModal();
    bottomPanel.style.display = 'none';

    // Atualiza UI e Gráfico
    updateCurrentProjectDisplay();
    updateApp();
    renderBottomPanel();

    Plotly.react(plotDiv, [], getPlotLayout(), getPlotConfig());
    fetchProjectsInfo();
}

// Carregar projeto salvo por nome
async function loadProjectByName(filename) {
    try {
        const res = await fetch(`/api/projects/load/${encodeURIComponent(filename)}`);
        if (!res.ok) throw new Error('Não foi possível ler o arquivo do servidor.');
        const resData = await res.json();
        if (resData.status === 'OK' && resData.data) {
            currentOpenProjectFilename = filename;
            currentOpenProjectPath = resData.filepath || null;
            await applyLoadedProjectData(resData.data);
            updateCurrentProjectDisplay();
            fetchProjectsInfo();
        }
    } catch (err) {
        alert('Erro ao carregar projeto: ' + err.message);
    }
}

// Aplica os dados do projeto no canvas
async function applyLoadedProjectData(data) {
    if (!data.samples || !Array.isArray(data.samples)) {
        alert('Formato de projeto inválido.');
        return;
    }

    useFidDomain();
    samples = data.samples.map((s, idx) => ({
        id: s.id || `sample_${idx}_${Date.now()}`,
        name: s.name || `Amostra ${idx + 1}`,
        color: s.color || COLOR_PALETTE[idx % COLOR_PALETTE.length],
        visible: s.visible !== undefined ? s.visible : true,
        n_points: s.n_points || s.time.length,
        sampling_interval_sec: s.sampling_interval_sec || 0.02,
        delay_time_sec: s.delay_time_sec || 0.0,
        shift: s.shift || 0,
        group: s.group || undefined,
        normManual: s.normManual > 0 ? s.normManual : 1,
        time: s.time || [],
        intensity: s.intensity || []
    }));
    compounds = Array.isArray(data.compounds) ? data.compounds.filter(c => c && c.id && c.name).map(c => ({
        id: String(c.id), name: String(c.name), rt: Number.isFinite(c.rt) ? c.rt : null,
        win: c.win > 0 ? c.win : 0.05,
        baseline: BASELINE_TYPES.includes(c.baseline) ? c.baseline : 'Valley-to-Valley',
        ion: c.ion > 0 ? c.ion : null, tol: c.tol > 0 ? c.tol : 0.5, is: !!c.is,
        domain: c.domain === 'ms' ? 'ms' : 'fid'
    })) : [];
    const ns = data.normSettings || {};
    normSettings = {
        mode: ['none', 'is', 'manual'].includes(ns.mode) ? ns.mode : 'none',
        isName: ns.isName || '',
        metric: ns.metric === 'height' ? 'height' : 'area',
        plot: !!ns.plot
    };

    appliedPeaks = data.appliedPeaks || [];
    resetPeakState();

    // GC-MS
    msSamples = Array.isArray(data.msSamples) ? data.msSamples.map((m, i) => ({
        ...m,
        color: m.color || COLOR_PALETTE[i % COLOR_PALETTE.length],
        scanIdx: Number.isInteger(m.scanIdx) ? m.scanIdx : Math.floor((m.time || []).length / 2),
        shift: m.shift || 0,
        group: m.group || undefined,
        normManual: m.normManual > 0 ? m.normManual : 1,
        avg: sanitizeMsAvg(m)
    })) : [];
    if (msSamples.some(m => m.shift)) {   // o servidor reinicia sem saber dos deslocamentos: reenvia
        fetch('/api/gcms/shifts', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ shifts: Object.fromEntries(msSamples.map(m => [m.id, m.shift || 0])) })
        }).catch(() => { });
    }
    msPanelMode = data.msPanelMode === 'avg' ? 'avg' : 'scan';
    msAvgResult = null;
    sampleGroups = Array.isArray(data.sampleGroups)
        ? data.sampleGroups.filter(g => g && g.id && g.name).map(g => ({ id: String(g.id), name: String(g.name), collapsed: !!g.collapsed, domain: g.domain === 'ms' ? 'ms' : 'fid' }))
        : [];
    pruneGroups();
    const savedEic = data.msEic;
    msEic = (savedEic && Array.isArray(savedEic.ions)) ? {
        ions: savedEic.ions.filter(x => Number.isFinite(x) && x > 0).slice(0, EIC_MAX),
        tol: (savedEic.tol > 0 && savedEic.tol <= 5) ? savedEic.tol : 0.5,
        showTic: !!savedEic.showTic
    } : { ions: [], tol: 0.5, showTic: false };
    msEicData = null;
    syncEicInputs();
    const mr = data.msReference;
    msReference = (mr && Array.isArray(mr.mz) && Array.isArray(mr.intensity) && mr.mz.length === mr.intensity.length && mr.mz.length > 0)
        ? { name: String(mr.name || 'Referência'), mz: mr.mz.map(Number), intensity: mr.intensity.map(Number), source: String(mr.source || 'projeto'), meta: mr.meta || {} }
        : null;
    msCmpResult = null;
    msCmpSource = 'scan';
    if (msPanelMode === 'cmp' && !msReference) msPanelMode = 'scan';
    activeMsId = msSamples.some(m => m.id === data.activeMsId) ? data.activeMsId : (msSamples[0] ? msSamples[0].id : null);
    msSpectrumCache.clear();
    msCurrentSpectrum = null;

    // Carrega Bloco de Notas e Rampas de Temperatura
    projectNotes = data.notes || '';
    if (projectNotesInput) projectNotesInput.value = projectNotes;

    if (Array.isArray(data.temperatureRamps) && data.temperatureRamps.length > 0) {
        temperatureRamps = data.temperatureRamps;
    } else {
        temperatureRamps = [
            {
                id: 'ramp_1',
                name: 'Rampa 1',
                visible: false,
                stages: [
                    { rate: null, temp: 40, hold: 2 },
                    { rate: 10, temp: 200, hold: 5 }
                ]
            }
        ];
    }
    selectedRampId = temperatureRamps[0]?.id || 'ramp_1';
    renderRampModal();

    if (data.currentTheme) {
        currentTheme = data.currentTheme;
        document.documentElement.setAttribute('data-theme', currentTheme);
        themeIcon.innerText = currentTheme === 'light' ? '🌙' : '☀️';
    }

    isYPercentMode = Boolean(data.isYPercentMode);
    currentOffsetY = Number(data.currentOffsetY) || 0;

    if (yModeRaw && yModePercent) {
        yModeRaw.checked = !isYPercentMode;
        yModePercent.checked = isYPercentMode;
    }
    if (ctxOffsetY) {
        ctxOffsetY.value = currentOffsetY;
    }

    await fetch('/api/project/sync', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(samples)
    });

    setMode(data.mode === 'gcms' ? 'gcms' : 'gcfid', true);

    Plotly.relayout(plotDiv, {
        'xaxis.autorange': true,
        'yaxis.autorange': !isYPercentMode,
        'yaxis.range': isYPercentMode ? [0, 105] : undefined,
        'yaxis.title': getYTitle()
    });

    Plotly.Plots.resize(plotDiv);
    setTimeout(() => Plotly.Plots.resize(plotDiv), 100);
    setTimeout(() => Plotly.Plots.resize(plotDiv), 300);
}

// ── MODAL DE INTEGRAÇÃO FORÇADA ──
function openIntegrationModal(initStart, initEnd) {
    setIntegrationModalMode(false);
    closeRemoveRegionModal(true);
    const ref = samples.find(s => s.visible) || samples[0];
    const maxPts = (ref.n_points || ref.time.length) - 1;

    [startSlider, endSlider, glStartSlider, glEndSlider, anchorSlider].forEach(slider => {
        slider.max = maxPts;
    });

    startSlider.value = initStart !== undefined ? initStart : Math.floor(maxPts * 0.15);
    endSlider.value = initEnd !== undefined ? initEnd : Math.floor(maxPts * 0.25);
    glStartSlider.value = Math.max(0, parseInt(startSlider.value) - 25);
    glEndSlider.value = Math.min(maxPts, parseInt(endSlider.value) + 25);
    anchorSlider.value = Math.min(maxPts, parseInt(endSlider.value) + 15);

    modal.style.display = 'flex';
    modalBodyContent.style.display = 'flex';
    updateSliderFeedbacks();
    renderCustomChecklist();
    updateApp();
}

function closeIntegrationModal() {
    modal.style.display = 'none';
    setIntegrationModalMode(false);
    updateApp();
}

// Alterna a janela entre "Integração Forçada" (nova) e "Editar integração" (pico existente)
function setIntegrationModalMode(isEdit, titleText) {
    const titleEl = modal.querySelector('.modal-title');
    const applySection = modal.querySelector('.apply-section');
    titleEl.innerText = isEdit ? titleText : '⏰ Integração Forçada';
    btnApplyIntegration.innerText = isEdit ? 'Salvar alteração' : 'Aplicar';
    applySection.style.display = isEdit ? 'none' : '';
    if (!isEdit) {
        editingPeak = null;
        document.getElementById('overlap-warning').style.display = 'none';
        document.getElementById('dropline-info').style.display = 'none';
    }
}

// Limites de um pico salvo, derivados da geometria guardada (funciona também em projetos antigos).
function derivePeakParams(pk, sample) {
    const maxPts = sample.time.length - 1;
    const m = pk.polygon.x.length / 2;   // o polígono é [x do pico] + [x do pico invertido]
    const start = findClosestIndex(sample, pk.polygon.x[0]);
    const end = findClosestIndex(sample, pk.polygon.x[m - 1]);
    const p = {
        start, end,
        glStart: Math.max(0, start - 25),
        glEnd: Math.min(maxPts, end + 25),
        anchor: Math.min(maxPts, end + 15)
    };
    if (pk.baseline_type === 'Drop Line') {
        p.glStart = findClosestIndex(sample, pk.baseline_line.x[0]);
        p.glEnd = findClosestIndex(sample, pk.baseline_line.x[1]);
    } else if (pk.baseline_type === 'Tangent Skim') {
        p.anchor = findClosestIndex(sample, pk.baseline_line.x[1]);
    }
    return p;
}

function openEditPeakModal(pk) {
    const sample = samples.find(s => s.id === pk.sample_id);
    if (!sample) {
        alert('A amostra deste pico não está mais no projeto.');
        return;
    }
    if (pickState !== 'idle') cancelPicking();
    closeRemoveRegionModal(true);
    if (!sample.visible) sample.visible = true;   // o preview só aparece com a amostra visível

    editingPeak = pk;
    const params = derivePeakParams(pk, sample);
    const maxPts = sample.time.length - 1;

    [startSlider, endSlider, glStartSlider, glEndSlider, anchorSlider].forEach(sl => { sl.max = maxPts; });
    startSlider.value = params.start;
    endSlider.value = params.end;
    glStartSlider.value = params.glStart;
    glEndSlider.value = params.glEnd;
    anchorSlider.value = params.anchor;

    baselineTypeSelect.value = pk.baseline_type;
    droplineControls.style.display = pk.baseline_type === 'Drop Line' ? 'block' : 'none';
    tangentControls.style.display = pk.baseline_type === 'Tangent Skim' ? 'block' : 'none';
    peakNameInput.value = pk.peak_name;

    const num = appliedPeaks.filter(p => p.sample_id === pk.sample_id).indexOf(pk) + 1;
    setIntegrationModalMode(true, `✏️ Editar integração — #${num} ${pk.peak_name} (${sample.name})`);

    modal.style.display = 'flex';
    modalBodyContent.style.display = 'flex';
    modalMinimizeBtn.innerText = '_';
    updateSliderFeedbacks();
}

function peakTimeRange(pk) {
    const m = pk.polygon.x.length / 2;
    return [pk.polygon.x[0], pk.polygon.x[m - 1]];
}

// Avisa (sem bloquear) quando os limites escolhidos invadem outro pico da mesma amostra
function updateOverlapWarning() {
    const box = document.getElementById('overlap-warning');
    const targets = getTargetSamples();
    const sIdx = parseInt(startSlider.value);
    const eIdx = parseInt(endSlider.value);
    const names = [];

    targets.forEach(s => {
        const n = s.time.length;
        const a = s.time[Math.max(0, Math.min(sIdx, n - 1))];
        const b = s.time[Math.max(0, Math.min(eIdx, n - 1))];
        if (!(a < b)) return;
        appliedPeaks.filter(p => p.sample_id === s.id).forEach((p, i) => {
            if (p === editingPeak) return;
            const [pa, pb] = peakTimeRange(p);
            if (Math.max(a, pa) < Math.min(b, pb)) {
                names.push(`#${i + 1} ${p.peak_name}` + (targets.length > 1 ? ` (${s.name})` : ''));
            }
        });
    });

    const lines = [];
    if (names.length > 0) {
        const shown = names.slice(0, 4).join(', ') + (names.length > 4 ? ` e mais ${names.length - 4}` : '');
        lines.push(`⚠ Os limites sobrepõem: ${shown}. A região em comum será contada nos dois picos.`);
    }
    // Drop Line em edição: nota discreta no ícone "i" ao lado do seletor (aparece ao passar o mouse)
    document.getElementById('dropline-info').style.display =
        (editingPeak && baselineTypeSelect.value === 'Drop Line') ? 'inline-flex' : 'none';

    box.innerHTML = '';
    lines.forEach(t => {
        const d = document.createElement('div');
        d.textContent = t;
        box.appendChild(d);
    });
    box.style.display = lines.length ? 'block' : 'none';
}

modalCloseBtn.addEventListener('click', closeIntegrationModal);
btnCancelIntegration.addEventListener('click', closeIntegrationModal);

modalMinimizeBtn.addEventListener('click', () => {
    const isMin = modalBodyContent.style.display === 'none';
    modalBodyContent.style.display = isMin ? 'flex' : 'none';
    modalMinimizeBtn.innerText = isMin ? '_' : '□';
});

// Modal Arrastável
let isDraggingModal = false, modalOffsetX = 0, modalOffsetY = 0;
modalDragHandle.addEventListener('mousedown', (e) => {
    if (e.target.closest('.modal-icon-btn')) return;
    isDraggingModal = true;
    modalOffsetX = e.clientX - modal.offsetLeft;
    modalOffsetY = e.clientY - modal.offsetTop;
});

document.addEventListener('mousemove', (e) => {
    if (!isDraggingModal) return;
    modal.style.left = `${e.clientX - modalOffsetX}px`;
    modal.style.top = `${e.clientY - modalOffsetY}px`;
    modal.style.right = 'auto';
});

document.addEventListener('mouseup', () => isDraggingModal = false);

// Função utilitária para arrastar modais
function makeDraggable(modalEl, handleEl) {
    if (!modalEl || !handleEl) return;
    let isDragging = false, offX = 0, offY = 0;
    handleEl.addEventListener('mousedown', (e) => {
        if (e.target.closest('.modal-icon-btn')) return;
        isDragging = true;
        offX = e.clientX - modalEl.offsetLeft;
        offY = e.clientY - modalEl.offsetTop;
    });
    document.addEventListener('mousemove', (e) => {
        if (!isDragging) return;
        modalEl.style.left = `${e.clientX - offX}px`;
        modalEl.style.top = `${e.clientY - offY}px`;
        modalEl.style.right = 'auto';
    });
    document.addEventListener('mouseup', () => {
        isDragging = false;
    });
}

// ── GERENCIAMENTO DO BLOCO DE NOTAS ──
function toggleNotesModal() {
    if (!notesModal) return;
    const isHidden = notesModal.style.display === 'none';
    notesModal.style.display = isHidden ? 'flex' : 'none';
    if (notesBodyContent) notesBodyContent.style.display = 'flex';
    if (notesBtn) notesBtn.classList.toggle('active-mode', isHidden);
    if (isHidden && projectNotesInput) {
        projectNotesInput.value = projectNotes || '';
        setTimeout(() => projectNotesInput.focus(), 60);
    }
}

if (notesBtn) notesBtn.addEventListener('click', toggleNotesModal);
if (notesCloseBtn) {
    notesCloseBtn.addEventListener('click', () => {
        notesModal.style.display = 'none';
        if (notesBtn) notesBtn.classList.remove('active-mode');
    });
}
if (notesMinimizeBtn) {
    notesMinimizeBtn.addEventListener('click', () => {
        const isMin = notesBodyContent.style.display === 'none';
        notesBodyContent.style.display = isMin ? 'flex' : 'none';
        notesMinimizeBtn.innerText = isMin ? '_' : '□';
    });
}
if (projectNotesInput) {
    projectNotesInput.addEventListener('input', () => {
        projectNotes = projectNotesInput.value;
    });
}
makeDraggable(notesModal, notesDragHandle);

// ── GERENCIAMENTO DAS RAMPAS DE TEMPERATURA ──
function toggleTempRampModal() {
    if (!tempRampModal) return;
    const isHidden = tempRampModal.style.display === 'none';
    tempRampModal.style.display = isHidden ? 'flex' : 'none';
    if (tempRampBodyContent) tempRampBodyContent.style.display = 'flex';
    if (tempRampBtn) tempRampBtn.classList.toggle('active-mode', isHidden);
    if (isHidden) {
        renderRampModal();
    }
}

if (tempRampBtn) tempRampBtn.addEventListener('click', toggleTempRampModal);
if (tempRampCloseBtn) {
    tempRampCloseBtn.addEventListener('click', () => {
        tempRampModal.style.display = 'none';
        if (tempRampBtn) tempRampBtn.classList.remove('active-mode');
    });
}
if (tempRampMinimizeBtn) {
    tempRampMinimizeBtn.addEventListener('click', () => {
        const isMin = tempRampBodyContent.style.display === 'none';
        tempRampBodyContent.style.display = isMin ? 'flex' : 'none';
        tempRampMinimizeBtn.innerText = isMin ? '_' : '□';
    });
}
makeDraggable(tempRampModal, tempRampDragHandle);

function getCurrentRamp() {
    if (!temperatureRamps || temperatureRamps.length === 0) {
        temperatureRamps = [
            {
                id: 'ramp_1',
                name: 'Rampa 1',
                visible: true,
                stages: [
                    { rate: null, temp: 40, hold: 2 },
                    { rate: 10, temp: 200, hold: 5 }
                ]
            }
        ];
        selectedRampId = temperatureRamps[0].id;
    }
    let ramp = temperatureRamps.find(r => r.id === selectedRampId);
    if (!ramp) {
        ramp = temperatureRamps[0];
        selectedRampId = ramp.id;
    }
    return ramp;
}

function updateRampSummary(ramp) {
    if (!rampSummaryInfo) return;
    const pts = calculateRampPoints(ramp, 0);
    const totalMin = pts.finalTime || 0;
    const initialTemp = ramp.stages && ramp.stages[0] ? ramp.stages[0].temp : 40;
    const finalTemp = pts.finalTemp !== undefined ? pts.finalTemp : initialTemp;
    rampSummaryInfo.innerText = `Tempo: ${totalMin.toFixed(1)} min | ${initialTemp} °C ➔ ${finalTemp} °C`;
}

function renderRampModal() {
    const ramp = getCurrentRamp();

    // Atualiza Select de Rampas
    if (rampSelect) {
        rampSelect.innerHTML = '';
        temperatureRamps.forEach(r => {
            const opt = document.createElement('option');
            opt.value = r.id;
            opt.innerText = `${r.name}${r.visible ? ' (Exibida)' : ''}`;
            if (r.id === selectedRampId) opt.selected = true;
            rampSelect.appendChild(opt);
        });
    }

    // Atualiza campos da rampa atual
    if (rampNameInput) rampNameInput.value = ramp.name || '';
    if (rampVisibleCheckbox) rampVisibleCheckbox.checked = Boolean(ramp.visible);

    // Renderiza tabela de estágios
    if (rampStagesBody) {
        rampStagesBody.innerHTML = '';
        const stages = ramp.stages || [];

        stages.forEach((st, idx) => {
            const tr = document.createElement('tr');
            if (idx === 0) {
                // Primeiro estágio: taxa é sempre '-'
                tr.innerHTML = `
                    <td style="font-weight:600; color:var(--text-muted);">${idx + 1}</td>
                    <td style="font-weight:bold; color:var(--text-muted); user-select:none;">—</td>
                    <td>
                        <input type="number" step="1" class="ramp-num-input" data-stage="0" data-field="temp" value="${st.temp !== undefined ? st.temp : 40}">
                    </td>
                    <td>
                        <input type="number" step="0.5" min="0" class="ramp-num-input" data-stage="0" data-field="hold" value="${st.hold !== undefined ? st.hold : 2}">
                    </td>
                    <td></td>
                `;
            } else {
                // Demais estágios: taxa configurável
                tr.innerHTML = `
                    <td style="font-weight:600; color:var(--text-muted);">${idx + 1}</td>
                    <td>
                        <input type="number" step="0.5" min="0.1" class="ramp-num-input" data-stage="${idx}" data-field="rate" value="${st.rate !== undefined && st.rate !== null ? st.rate : 10}">
                    </td>
                    <td>
                        <input type="number" step="1" class="ramp-num-input" data-stage="${idx}" data-field="temp" value="${st.temp !== undefined ? st.temp : 200}">
                    </td>
                    <td>
                        <input type="number" step="0.5" min="0" class="ramp-num-input" data-stage="${idx}" data-field="hold" value="${st.hold !== undefined ? st.hold : 5}">
                    </td>
                    <td>
                        <button class="btn-del-stage" data-stage="${idx}" title="Excluir estágio">✕</button>
                    </td>
                `;
            }
            rampStagesBody.appendChild(tr);
        });
    }

    updateRampSummary(ramp);
}

// Eventos de Rampa
if (rampSelect) {
    rampSelect.addEventListener('change', () => {
        selectedRampId = rampSelect.value;
        renderRampModal();
        updateApp();
    });
}

if (rampNameInput) {
    rampNameInput.addEventListener('input', () => {
        const ramp = getCurrentRamp();
        ramp.name = rampNameInput.value.trim() || 'Rampa';
        if (rampSelect) {
            const opt = rampSelect.querySelector(`option[value="${ramp.id}"]`);
            if (opt) opt.innerText = `${ramp.name}${ramp.visible ? ' (Exibida)' : ''}`;
        }
        updateApp();
    });
}

if (rampVisibleCheckbox) {
    rampVisibleCheckbox.addEventListener('change', () => {
        const ramp = getCurrentRamp();
        ramp.visible = rampVisibleCheckbox.checked;
        if (rampSelect) {
            const opt = rampSelect.querySelector(`option[value="${ramp.id}"]`);
            if (opt) opt.innerText = `${ramp.name}${ramp.visible ? ' (Exibida)' : ''}`;
        }
        updateApp();
    });
}

if (btnAddRamp) {
    btnAddRamp.addEventListener('click', () => {
        const newId = 'ramp_' + Date.now();
        const newRamp = {
            id: newId,
            name: `Rampa ${temperatureRamps.length + 1}`,
            visible: true,
            stages: [
                { rate: null, temp: 40, hold: 2 },
                { rate: 10, temp: 240, hold: 5 }
            ]
        };
        temperatureRamps.push(newRamp);
        selectedRampId = newId;
        renderRampModal();
        updateApp();
    });
}

if (btnDelRamp) {
    btnDelRamp.addEventListener('click', () => {
        const ramp = getCurrentRamp();
        if (temperatureRamps.length > 1) {
            if (confirm(`Deseja excluir "${ramp.name}"?`)) {
                temperatureRamps = temperatureRamps.filter(r => r.id !== ramp.id);
                selectedRampId = temperatureRamps[0].id;
                renderRampModal();
                updateApp();
            }
        } else {
            ramp.name = 'Rampa 1';
            ramp.visible = false;
            ramp.stages = [
                { rate: null, temp: 40, hold: 2 },
                { rate: 10, temp: 200, hold: 5 }
            ];
            renderRampModal();
            updateApp();
        }
    });
}

if (btnAddStage) {
    btnAddStage.addEventListener('click', () => {
        const ramp = getCurrentRamp();
        const lastSt = ramp.stages[ramp.stages.length - 1];
        const lastTemp = lastSt ? Number(lastSt.temp) || 200 : 200;
        ramp.stages.push({
            rate: 10,
            temp: Math.min(450, lastTemp + 40),
            hold: 2
        });
        renderRampModal();
        updateApp();
    });
}

if (rampStagesBody) {
    rampStagesBody.addEventListener('input', (e) => {
        const inp = e.target;
        if (!inp.classList.contains('ramp-num-input')) return;
        const stageIdx = parseInt(inp.dataset.stage);
        const field = inp.dataset.field;
        const val = parseFloat(inp.value);
        const ramp = getCurrentRamp();
        if (ramp && ramp.stages && ramp.stages[stageIdx]) {
            ramp.stages[stageIdx][field] = isNaN(val) ? 0 : val;
            updateRampSummary(ramp);
            updateApp();
        }
    });

    rampStagesBody.addEventListener('click', (e) => {
        const btn = e.target.closest('.btn-del-stage');
        if (!btn) return;
        const stageIdx = parseInt(btn.dataset.stage);
        const ramp = getCurrentRamp();
        if (ramp && stageIdx > 0 && ramp.stages.length > 1) {
            ramp.stages.splice(stageIdx, 1);
            renderRampModal();
            updateApp();
        }
    });
}

function updateSliderFeedbacks() {
    if (samples.length === 0) return;
    const ref = (editingPeak && samples.find(s => s.id === editingPeak.sample_id)) || samples[0];
    const dt = ref.sampling_interval_sec || 0.02;
    const total = ref.n_points || ref.time.length;

    const sVal = parseInt(startSlider.value);
    const eVal = parseInt(endSlider.value);

    startIdxFeedback.innerText = `${sVal} / ${total}`;
    startTimeFeedback.innerText = `${(sVal * dt).toFixed(4)} s`;

    endIdxFeedback.innerText = `${eVal} / ${total}`;
    endTimeFeedback.innerText = `${(eVal * dt).toFixed(4)} s`;

    const glSVal = parseInt(glStartSlider.value);
    const glEVal = parseInt(glEndSlider.value);
    document.getElementById('gl-start-idx-feedback').innerText = `${glSVal} / ${total}`;
    document.getElementById('gl-start-time-feedback').innerText = `${(glSVal * dt).toFixed(4)} s`;
    document.getElementById('gl-end-idx-feedback').innerText = `${glEVal} / ${total}`;
    document.getElementById('gl-end-time-feedback').innerText = `${(glEVal * dt).toFixed(4)} s`;

    const ancVal = parseInt(anchorSlider.value);
    document.getElementById('anchor-idx-feedback').innerText = `${ancVal} / ${total}`;
    document.getElementById('anchor-time-feedback').innerText = `${(ancVal * dt).toFixed(4)} s`;

    updateApp();
}

function bindNudgeButtons(prefix, slider) {
    const b10 = document.getElementById(`${prefix}-step-back-10`);
    const b1 = document.getElementById(`${prefix}-step-back-1`);
    const f1 = document.getElementById(`${prefix}-step-fwd-1`);
    const f10 = document.getElementById(`${prefix}-step-fwd-10`);

    if (b10) b10.addEventListener('click', () => { slider.value = Math.max(0, parseInt(slider.value) - 10); updateSliderFeedbacks(); });
    if (b1) b1.addEventListener('click', () => { slider.value = Math.max(0, parseInt(slider.value) - 1); updateSliderFeedbacks(); });
    if (f1) f1.addEventListener('click', () => { slider.value = Math.min(parseInt(slider.max), parseInt(slider.value) + 1); updateSliderFeedbacks(); });
    if (f10) f10.addEventListener('click', () => { slider.value = Math.min(parseInt(slider.max), parseInt(slider.value) + 10); updateSliderFeedbacks(); });

    slider.addEventListener('input', updateSliderFeedbacks);
}

bindNudgeButtons('start', startSlider);
bindNudgeButtons('end', endSlider);
bindNudgeButtons('gl-start', glStartSlider);
bindNudgeButtons('gl-end', glEndSlider);
bindNudgeButtons('anchor', anchorSlider);

baselineTypeSelect.addEventListener('change', () => {
    const val = baselineTypeSelect.value;
    droplineControls.style.display = val === 'Drop Line' ? 'block' : 'none';
    tangentControls.style.display = val === 'Tangent Skim' ? 'block' : 'none';
    updateApp();
});

targetVisible.addEventListener('change', () => { customSelectionBox.style.display = 'none'; updateApp(); });
targetCustom.addEventListener('change', () => { customSelectionBox.style.display = 'flex'; updateApp(); });

function renderCustomChecklist() {
    customChecklist.innerHTML = '';
    samples.forEach(s => {
        const li = document.createElement('li');
        li.className = 'check-item';
        li.innerHTML = `
            <input type="checkbox" value="${escapeHtml(s.id)}" checked>
            <span class="color-dot" style="background-color: ${s.color}; width: 8px; height: 8px;"></span>
            <span>${escapeHtml(s.name)}</span>
        `;
        li.querySelector('input').addEventListener('change', updateApp);
        customChecklist.appendChild(li);
    });
}

btnSelectAll.addEventListener('click', () => { customChecklist.querySelectorAll('input').forEach(cb => cb.checked = true); updateApp(); });
btnDeselectAll.addEventListener('click', () => { customChecklist.querySelectorAll('input').forEach(cb => cb.checked = false); updateApp(); });

// ── AVISO RÁPIDO (toast) ──
let toastTimer = null;
function showToast(msg) {
    let el = document.getElementById('app-toast');
    if (!el) {
        el = document.createElement('div');
        el.id = 'app-toast';
        // Estilos inline de propósito: o aviso não pode depender do style.css (o navegador pode ter uma versão antiga em cache)
        Object.assign(el.style, {
            position: 'fixed', left: '50%', bottom: '28px', transform: 'translateX(-50%)',
            zIndex: '10000', maxWidth: '80vw', padding: '9px 16px', borderRadius: '8px',
            fontSize: '13px', fontFamily: 'inherit', boxShadow: '0 6px 20px rgba(0,0,0,0.3)',
            pointerEvents: 'none', opacity: '0', transition: 'opacity 0.18s ease'
        });
        document.body.appendChild(el);
    }
    const dark = currentTheme === 'dark';
    el.style.background = dark ? '#e2e8f0' : '#0f172a';
    el.style.color = dark ? '#0f172a' : '#f8fafc';
    el.textContent = msg;
    void el.offsetWidth;   // garante a transição mesmo com avisos em sequência
    el.style.opacity = '1';
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => { el.style.opacity = '0'; }, 2400);
}

// ── DESFAZER / REFAZER (Ctrl+Z, Ctrl+Y ou Ctrl+Shift+Z) — somente atalhos, sem botões ──
// Cobre: integrar, editar integração, renomear e excluir picos.
// Limpar/fechar projeto, abrir projeto e remover amostra zeram o histórico.
function pushHistory(label) {
    undoStack.push({ label, peaks: appliedPeaks.slice() });
    if (undoStack.length > HISTORY_LIMIT) undoStack.shift();
    redoStack = [];
}

function resetPeakState() {   // histórico de desfazer/refazer + linhas marcadas na tabela
    undoStack = [];
    redoStack = [];
    selectedPeaks.clear();
}

function restorePeaks(peaks) {
    appliedPeaks = peaks.slice();
    selectedPeaks.forEach(p => { if (!appliedPeaks.includes(p)) selectedPeaks.delete(p); });
    // se a aba ativa ficou sem picos, volta para "Todas"
    if (activePanelTab !== 'all' && !appliedPeaks.some(p => p.sample_id === activePanelTab)) activePanelTab = 'all';
    updateApp();
    renderBottomPanel();
}

function undoPeaks() {
    if (undoStack.length === 0) { showToast('Nada para desfazer'); return; }
    const entry = undoStack.pop();
    redoStack.push({ label: entry.label, peaks: appliedPeaks.slice() });
    restorePeaks(entry.peaks);
    showToast(`Desfeito: ${entry.label}`);
}

function redoPeaks() {
    if (redoStack.length === 0) { showToast('Nada para refazer'); return; }
    const entry = redoStack.pop();
    undoStack.push({ label: entry.label, peaks: appliedPeaks.slice() });
    restorePeaks(entry.peaks);
    showToast(`Refeito: ${entry.label}`);
}

document.addEventListener('keydown', (e) => {
    if (!(e.ctrlKey || e.metaKey) || e.altKey) return;
    const key = e.code === 'KeyZ' ? 'z' : e.code === 'KeyY' ? 'y' : e.key.toLowerCase();
    if (key !== 'z' && key !== 'y') return;

    // Em campos de texto, o desfazer nativo do navegador continua valendo
    const tag = (e.target.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'textarea' || tag === 'select' || e.target.isContentEditable) return;

    if (modal.style.display !== 'none' || removeModal.style.display !== 'none' || pickState !== 'idle') return;   // ação em andamento

    e.preventDefault();
    if (key === 'y' || (key === 'z' && e.shiftKey)) redoPeaks(); else undoPeaks();
});

// ── APLICAR INTEGRAÇÃO VIA BACKEND ──
btnApplyIntegration.addEventListener('click', async () => {
    const editing = editingPeak;
    const targetSamples = getTargetSamples();
    const targetIds = targetSamples.map(s => s.id);

    if (targetIds.length === 0) {
        alert('Nenhuma amostra selecionada.');
        return;
    }

    const purePeakName = peakNameInput.value.trim() || (editing ? editing.peak_name : `Pico ${appliedPeaks.length + 1}`);

    const payload = {
        sample_ids: targetIds,
        peak_name: purePeakName,
        baseline_type: baselineTypeSelect.value,
        start_idx: parseInt(startSlider.value),
        end_idx: parseInt(endSlider.value),
        gl_start_idx: parseInt(glStartSlider.value),
        gl_end_idx: parseInt(glEndSlider.value),
        anchor_idx: parseInt(anchorSlider.value)
    };

    btnApplyIntegration.innerText = 'Calculando...';
    btnApplyIntegration.disabled = true;

    try {
        const res = await fetch('/api/integrate/forced', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });

        const data = await res.json();
        if (data.status === 'OK' && data.peaks) {
            if (editing) {
                const idx = appliedPeaks.indexOf(editing);
                if (idx === -1 || data.peaks.length === 0) {
                    alert('Não foi possível recalcular: verifique os limites ou se o pico ainda existe.');
                    return;
                }
                pushHistory(`Editar integração de "${editing.peak_name}"`);
                replacePeakAt(idx, data.peaks[0]);   // mantém a posição (e a numeração) do pico
            } else {
                pushHistory(`Integrar "${purePeakName}"` + (data.peaks.length > 1 ? ` (${data.peaks.length} amostras)` : ''));
                appliedPeaks = appliedPeaks.concat(data.peaks);
            }
            closeIntegrationModal();
            peakNameInput.value = '';
            bottomPanel.style.display = 'flex';
            renderBottomPanel();
        }
    } catch (err) {
        alert('Erro ao calcular no backend.');
    } finally {
        btnApplyIntegration.innerText = editingPeak ? 'Salvar alteração' : 'Aplicar';
        btnApplyIntegration.disabled = false;
    }
});

// ── EXPORTAR TABELA EM CSV ──
btnExportTableCsv.addEventListener('click', () => {
    let visiblePeaks = getVisiblePeaks();
    let tabName = "Todas";

    if (activePanelTab !== 'all') {
        const s = samples.find(item => item.id === activePanelTab);
        tabName = s ? s.name : activePanelTab;
    }

    if (visiblePeaks.length === 0) {
        alert('Não há dados nesta guia para exportar.');
        return;
    }

    let csvContent = "\uFEFF";
    const normOn = normSettings.mode !== 'none';
    csvContent += "#;Amostra;Nome;tR (min);Area (u.a.*s);Altura;% Area;Baseline" + (normOn ? `;${normSettings.metric === 'height' ? 'Altura' : 'Area'} norm.` : '') + "\n";

    visiblePeaks.forEach(pk => {
        const parent = samples.find(s => s.id === pk.sample_id);
        const runPeaks = appliedPeaks.filter(p => p.sample_id === pk.sample_id);
        const peakNum = runPeaks.indexOf(pk) + 1;
        const totalAreaRun = getAreaPctBase(pk);
        const areaPct = totalAreaRun > 0 ? ((pk.area / totalAreaRun) * 100).toFixed(6) : "100.000000";

        const row = [
            peakNum,
            parent ? parent.name : pk.sample_id,
            `"${pk.peak_name}"`,
            pk.rt_min.toString(),
            pk.area.toString(),
            pk.height.toString(),
            showAreaPctColumn() ? areaPct : "",
            pk.baseline_type
        ];
        if (normOn) row.push(getNormalizedValue(pk) === null ? '' : String(getNormalizedValue(pk)));
        csvContent += row.join(";") + "\n";
    });

    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `integracao_${tabName.replace(/\s+/g, '_')}_${new Date().toISOString().slice(0, 10)}.csv`;
    a.click();
    URL.revokeObjectURL(url);
});

// ── GAVETA INFERIOR DE RESULTADOS ──
function renderBottomPanel() {
    const peaksHere = domainPeaks();
    if (peaksHere.length === 0) {
        bottomPanel.style.display = 'none';
        return;
    }
    if (activePanelTab !== 'all' && !peaksHere.some(p => p.sample_id === activePanelTab)) activePanelTab = 'all';
    bottomPanel.style.display = 'flex';

    panelTabList.innerHTML = '';

    const allTab = document.createElement('button');
    allTab.className = `tab-btn ${activePanelTab === 'all' ? 'active' : ''}`;
    allTab.innerText = `Todas (${peaksHere.length})`;
    allTab.addEventListener('click', () => { activePanelTab = 'all'; selectedPeaks.clear(); renderBottomPanel(); });
    panelTabList.appendChild(allTab);

    const sampleIdsWithPeaks = [...new Set(peaksHere.map(p => p.sample_id))];
    sampleIdsWithPeaks.forEach(sId => {
        const sample = samples.find(s => s.id === sId);
        if (!sample) return;
        const count = appliedPeaks.filter(p => p.sample_id === sId).length;

        const tab = document.createElement('button');
        tab.className = `tab-btn ${activePanelTab === sId ? 'active' : ''}`;
        tab.innerHTML = `<span style="color:${sample.color}">●</span> ${escapeHtml(sample.name)} (${count})`;
        tab.addEventListener('click', () => { activePanelTab = sId; selectedPeaks.clear(); renderBottomPanel(); });
        panelTabList.appendChild(tab);
    });

    thAreaPct.style.display = showAreaPctColumn() ? 'table-cell' : 'none';
    const normOn = normSettings.mode !== 'none';
    thNorm.style.display = normOn ? 'table-cell' : 'none';
    thNorm.textContent = normSettings.metric === 'height' ? 'Altura norm.' : 'Área norm.';

    let visiblePeaks = peaksHere;
    if (activePanelTab !== 'all') {
        visiblePeaks = peaksHere.filter(p => p.sample_id === activePanelTab);
    }

    peaksTableBody.innerHTML = '';
    visiblePeaks.forEach((pk) => {
        const parent = samples.find(s => s.id === pk.sample_id);
        const runPeaks = appliedPeaks.filter(p => p.sample_id === pk.sample_id);
        const peakNum = runPeaks.indexOf(pk) + 1;

        const totalAreaRun = getAreaPctBase(pk);
        const areaPct = totalAreaRun > 0 ? ((pk.area / totalAreaRun) * 100).toFixed(2) + ' %' : '100.00 %';

        const tr = document.createElement('tr');
        tr.innerHTML = `
            <td><input type="checkbox" class="peak-select" ${selectedPeaks.has(pk) ? 'checked' : ''} title="Marcar para remover"></td>
            <td><strong>${peakNum}</strong></td>
            <td><span style="color:${parent ? parent.color : '#fff'}">●</span> ${escapeHtml(parent ? parent.name : pk.sample_id)}</td>
            <td><span class="editable-name" title="Duplo clique para renomear">${escapeHtml(pk.peak_name)}</span></td>
            <td>${pk.rt_min.toFixed(4)}</td>
            <td title="Valor completo: ${pk.area}">
                ${pk.area.toLocaleString(undefined, { minimumFractionDigits: 2, maximumFractionDigits: 2 })}
                <button class="btn-copy" title="Copiar área com todas as casas decimais">📋</button>
            </td>
            <td title="Valor completo: ${pk.height}">${pk.height.toLocaleString(undefined, { minimumFractionDigits: 2, maximumFractionDigits: 2 })}</td>
            <td style="display: ${normOn ? 'table-cell' : 'none'}" title="Normalizado conforme o modo escolhido em ⇄">${formatNormValue(pk)}</td>
            <td style="display: ${showAreaPctColumn() ? 'table-cell' : 'none'}"${appMode === 'gcms' && activePanelTab === 'all' ? ' title="Porcentagem da soma das áreas de todos os picos exibidos desta corrida"' : ''}><strong>${areaPct}</strong></td>
            <td><small class="baseline-edit" style="color:var(--text-muted)" title="Clique para editar a integração deste pico">${escapeHtml(pk.baseline_type)} ✎</small></td>
            <td><button class="btn-close" title="Excluir este pico">✕</button></td>
        `;

        tr.querySelector('.btn-copy').addEventListener('click', (e) => {
            e.stopPropagation();
            navigator.clipboard.writeText(pk.area.toString());
            const btn = e.currentTarget;
            btn.innerText = '✓';
            setTimeout(() => btn.innerText = '📋', 1200);
        });

        tr.querySelector('.baseline-edit').addEventListener('click', () => openEditPeakModal(pk));

        const rowCheck = tr.querySelector('.peak-select');
        tr.classList.toggle('selected-row', rowCheck.checked);
        rowCheck.addEventListener('change', () => {
            if (rowCheck.checked) selectedPeaks.add(pk); else selectedPeaks.delete(pk);
            tr.classList.toggle('selected-row', rowCheck.checked);
            updateTableSelectionUI();
        });

        const nameEl = tr.querySelector('.editable-name');
        nameEl.addEventListener('dblclick', () => {
            const currentName = pk.peak_name;
            const input = document.createElement('input');
            input.type = 'text';
            input.className = 'rename-input';
            input.value = currentName;
            nameEl.replaceWith(input);
            input.focus();
            input.select();

            let cancelled = false;
            const finishRename = () => {
                const newName = input.value.trim() || currentName;
                const idx = appliedPeaks.indexOf(pk);
                if (!cancelled && idx > -1 && newName !== currentName) {
                    pushHistory(`Renomear pico "${currentName}"`);
                    replacePeakAt(idx, { ...pk, peak_name: newName });   // novo objeto (histórico)
                }
                updateApp();
                renderBottomPanel();
            };

            input.addEventListener('keydown', (ke) => {
                if (ke.key === 'Enter') finishRename();
                if (ke.key === 'Escape') { cancelled = true; renderBottomPanel(); }
            });
            input.addEventListener('blur', finishRename);
        });

        tr.querySelector('.btn-close').addEventListener('click', () => {
            removePeaks(new Set([pk]), `Excluir pico "${pk.peak_name}"`);
        });

        peaksTableBody.appendChild(tr);
    });

    updateTableSelectionUI();
}

btnToggleShading.addEventListener('click', () => {
    isShadingVisible = !isShadingVisible;
    btnToggleShading.classList.toggle('off', !isShadingVisible);
    updateApp();
});

btnToggleLabels.addEventListener('click', () => {
    isLabelsVisible = !isLabelsVisible;
    btnToggleLabels.classList.toggle('off', !isLabelsVisible);
    updateApp();
});

// Redimensionamento Vertical da Gaveta
let isResizingBottom = false;
panelResizer.addEventListener('mousedown', () => {
    isResizingBottom = true;
    panelResizer.classList.add('resizing');
    document.body.style.cursor = 'row-resize';
    document.body.style.userSelect = 'none';
});

document.addEventListener('mousemove', (e) => {
    if (!isResizingBottom) return;
    const newHeight = window.innerHeight - e.clientY;
    if (newHeight >= 80 && newHeight <= window.innerHeight * 0.75) {
        bottomPanel.classList.remove('minimized');
        bottomPanel.style.height = `${newHeight}px`;
        Plotly.Plots.resize(plotDiv);
    }
});

document.addEventListener('mouseup', () => {
    if (isResizingBottom) {
        isResizingBottom = false;
        panelResizer.classList.remove('resizing');
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
        Plotly.Plots.resize(plotDiv);
    }
});

toggleTableBtn.addEventListener('click', () => {
    if (bottomPanel.style.display === 'none') {
        bottomPanel.style.display = 'flex';
        renderBottomPanel();
    } else {
        bottomPanel.style.display = 'none';
    }
    Plotly.Plots.resize(plotDiv);
});

panelMinimizeBtn.addEventListener('click', () => {
    const isMin = bottomPanel.classList.toggle('minimized');
    panelMinimizeBtn.innerText = isMin ? '□' : '_';
    Plotly.Plots.resize(plotDiv);
});

panelCloseBtn.addEventListener('click', () => {
    bottomPanel.style.display = 'none';
    Plotly.Plots.resize(plotDiv);
});

// ── COLAPSAR E EXPANDIR BARRA LATERAL ESQUERDA ──
collapseBtn.addEventListener('click', () => {
    sidebar.classList.add('collapsed');
    resizer.classList.add('hidden');
    expandBtn.style.display = 'flex';
    setTimeout(() => Plotly.Plots.resize(plotDiv), 180);
});

expandBtn.addEventListener('click', () => {
    sidebar.classList.remove('collapsed');
    resizer.classList.remove('hidden');
    expandBtn.style.display = 'none';
    setTimeout(() => Plotly.Plots.resize(plotDiv), 180);
});

// ── UPLOAD DE ARQUIVOS .CDF ──
uploadBtn.addEventListener('click', () => fileInput.click());
fileInput.addEventListener('change', (e) => {
    if (e.target.files.length > 0) uploadFiles(e.target.files);
});

sidebar.addEventListener('dragover', (e) => { e.preventDefault(); sidebar.classList.add('dragover'); });
sidebar.addEventListener('dragleave', (e) => { if (!sidebar.contains(e.relatedTarget)) sidebar.classList.remove('dragover'); });
sidebar.addEventListener('drop', (e) => {
    e.preventDefault();
    sidebar.classList.remove('dragover');
    if (e.dataTransfer.files.length > 0) uploadFiles(e.dataTransfer.files);
});

async function uploadFiles(files) {
    if (appMode === 'gcms') return uploadMsFiles(files);
    const formData = new FormData();
    for (let i = 0; i < files.length; i++) formData.append('files', files[i]);

    uploadBtn.disabled = true;
    try {
        const res = await fetch('/api/upload/gc-fid', { method: 'POST', body: formData });
        const data = await res.json();
        if (data.success && data.success.length > 0) {
            data.success.forEach(chrom => {
                const color = COLOR_PALETTE[samples.length % COLOR_PALETTE.length];
                samples.push({
                    id: chrom.sample_id,
                    name: sanitizeFilename(chrom.filename),
                    n_points: chrom.n_points,
                    sampling_interval_sec: chrom.sampling_interval_sec,
                    delay_time_sec: chrom.delay_time_sec || 0.0,
                    time: chrom.time,
                    intensity: chrom.intensity,
                    color: color,
                    visible: true
                });
            });
            updateApp();
        }
    } catch (err) {
        alert('Erro ao carregar arquivos.');
    } finally {
        uploadBtn.disabled = false;
        fileInput.value = '';
    }
}

clearBtn.addEventListener('click', () => {
    useFidDomain();
    samples = [];
    appliedPeaks = [];
    resetPeakState();
    clearMsState();
    activateDomain();
    normSettings = { mode: 'none', isName: '', metric: 'area', plot: false };
    compounds = [];
    sampleGroups = [];
    currentOpenProjectFilename = null;
    currentOpenProjectPath = null;
    projectNotes = '';
    if (projectNotesInput) projectNotesInput.value = '';
    temperatureRamps = [
        {
            id: 'ramp_1',
            name: 'Rampa 1',
            visible: false,
            stages: [
                { rate: null, temp: 40, hold: 2 },
                { rate: 10, temp: 200, hold: 5 }
            ]
        }
    ];
    selectedRampId = 'ramp_1';
    renderRampModal();
    if (notesModal) notesModal.style.display = 'none';
    if (tempRampModal) tempRampModal.style.display = 'none';
    if (notesBtn) notesBtn.classList.remove('active-mode');
    if (tempRampBtn) tempRampBtn.classList.remove('active-mode');
    cancelPicking();
    closeIntegrationModal();
    bottomPanel.style.display = 'none';
    updateCurrentProjectDisplay();
    updateApp();
    fetchProjectsInfo();
});

themeBtn.addEventListener('click', () => {
    currentTheme = currentTheme === 'light' ? 'dark' : 'light';
    document.documentElement.setAttribute('data-theme', currentTheme);
    themeIcon.innerText = currentTheme === 'light' ? '🌙' : '☀️';
    if (appMode === 'gcms') {
        updateApp();
        renderSpectrum();
    } else {
        Plotly.relayout(plotDiv, getPlotLayout());
    }
});

// Redimensionamento Barra Esquerda
let isResizing = false;
resizer.addEventListener('mousedown', () => {
    isResizing = true;
    resizer.classList.add('resizing');
    document.body.style.cursor = 'col-resize';
    document.body.style.userSelect = 'none';
});

document.addEventListener('mousemove', (e) => {
    if (!isResizing) return;
    const newWidth = e.clientX;
    if (newWidth >= 180 && newWidth <= 500) {
        sidebar.style.width = `${newWidth}px`;
        Plotly.Plots.resize(plotDiv);
    }
});

document.addEventListener('mouseup', () => {
    if (isResizing) {
        isResizing = false;
        resizer.classList.remove('resizing');
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
        Plotly.Plots.resize(plotDiv);
    }
});

// Redimensionamento Barra Direita (Projetos)
let isResizingRight = false;
rightResizer.addEventListener('mousedown', () => {
    isResizingRight = true;
    rightResizer.classList.add('resizing');
    document.body.style.cursor = 'col-resize';
    document.body.style.userSelect = 'none';
});

document.addEventListener('mousemove', (e) => {
    if (!isResizingRight) return;
    const newWidth = window.innerWidth - e.clientX;
    if (newWidth >= 180 && newWidth <= 500) {
        rightSidebar.style.width = `${newWidth}px`;
        Plotly.Plots.resize(plotDiv);
    }
});

document.addEventListener('mouseup', () => {
    if (isResizingRight) {
        isResizingRight = false;
        rightResizer.classList.remove('resizing');
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
        Plotly.Plots.resize(plotDiv);
    }
});

// ══════════════════════════════════════════════════════════════
// ── REMOÇÃO DE PICOS ──
// Dois caminhos, uma única rotina (removePeaks), sempre desfazível com Ctrl+Z:
//   1) botão 🧹 da barra esquerda: escolhe início/fim no gráfico e marca as corridas na janela;
//   2) botão da gaveta (tabela): remove os picos MARCADOS; sem marcação, todos os da guia ativa.
// ══════════════════════════════════════════════════════════════
const regionBtn = document.getElementById('region-remove-btn');
const removeModal = document.getElementById('remove-region-modal');
const regionT0 = document.getElementById('region-t0');
const regionT1 = document.getElementById('region-t1');
const regionRunList = document.getElementById('region-run-list');
const btnConfirmRegionRemove = document.getElementById('btn-confirm-region-remove');
const btnRemovePeaks = document.getElementById('btn-remove-peaks');
const peaksSelectAll = document.getElementById('peaks-select-all');
let regionUnchecked = new Set();   // corridas que o usuário desmarcou na janela (o padrão é marcar todas com picos)

const plural = (n, word) => `${n} ${word}${n === 1 ? '' : 's'}`;

function replacePeakAt(idx, newPk) {
    const old = appliedPeaks[idx];
    appliedPeaks[idx] = newPk;
    if (selectedPeaks.delete(old)) selectedPeaks.add(newPk);   // a marcação acompanha o pico editado/renomeado
}

function removePeaks(set, label) {
    if (set.size === 0) return;
    pushHistory(label);
    appliedPeaks = appliedPeaks.filter(p => !set.has(p));
    set.forEach(p => selectedPeaks.delete(p));
    if (activePanelTab !== 'all' && !appliedPeaks.some(p => p.sample_id === activePanelTab)) activePanelTab = 'all';
    updateApp();
    renderBottomPanel();
    showToast(`${plural(set.size, 'pico')} removido${set.size === 1 ? '' : 's'} — Ctrl+Z desfaz`);
}

// ── 1) Por região ──
regionBtn.addEventListener('click', toggleRegionRemoveMode);
makeDraggable(removeModal, document.getElementById('remove-region-drag-handle'));

function getRegionBounds() {
    const a = parseFloat(regionT0.value), b = parseFloat(regionT1.value);
    if (isNaN(a) || isNaN(b)) return null;
    return a <= b ? [a, b] : [b, a];
}

function getRegionShapes() {
    if (!regionRemoval || regionRemoval.t0 === null) return [];
    return [{
        type: 'rect', xref: 'x', yref: 'paper',
        x0: regionRemoval.t0, x1: regionRemoval.t1, y0: 0, y1: 1,
        fillcolor: 'rgba(239, 68, 68, 0.10)',
        line: { color: 'rgba(239, 68, 68, 0.65)', width: 1, dash: 'dot' },
        layer: 'below'
    }];
}

function openRemoveRegionModal(t0, t1) {
    if (modal.style.display !== 'none') { modal.style.display = 'none'; setIntegrationModalMode(false); }
    regionUnchecked = new Set();
    regionT0.value = t0.toFixed(3);
    regionT1.value = t1.toFixed(3);
    removeModal.style.display = 'flex';
    refreshRegionModal();
}

// Recalcula, para cada corrida, os picos cujo ápice cai na região (e o destaque no gráfico)
function refreshRegionModal() {
    const bounds = getRegionBounds();
    const doomed = new Set();
    regionRunList.innerHTML = '';

    const withPeaks = new Set(appliedPeaks.map(p => p.sample_id));
    samples.filter(s => withPeaks.has(s.id)).forEach(s => {
        const inRegion = bounds
            ? appliedPeaks.filter(p => p.sample_id === s.id && p.rt_min >= bounds[0] && p.rt_min <= bounds[1])
            : [];
        const enabled = inRegion.length > 0;
        const checked = enabled && !regionUnchecked.has(s.id);
        if (checked) inRegion.forEach(p => doomed.add(p));

        const li = document.createElement('li');
        li.className = 'check-item region-run' + (enabled ? '' : ' disabled');
        li.innerHTML = `
            <input type="checkbox" ${checked ? 'checked' : ''} ${enabled ? '' : 'disabled'}>
            <span class="color-dot" style="background-color: ${s.color}; width: 8px; height: 8px;"></span>
            <span class="region-run-name">${escapeHtml(s.name)}</span>
            <span class="region-run-count" title="${escapeHtml(inRegion.map(p => p.peak_name).join(', '))}">${enabled ? plural(inRegion.length, 'pico') : 'sem picos'}</span>
        `;
        li.querySelector('input').addEventListener('change', (e) => {
            if (e.target.checked) regionUnchecked.delete(s.id); else regionUnchecked.add(s.id);
            refreshRegionModal();
        });
        regionRunList.appendChild(li);
    });

    regionRemoval = { t0: bounds ? bounds[0] : null, t1: bounds ? bounds[1] : null, doomed };
    btnConfirmRegionRemove.textContent = doomed.size ? `Remover ${plural(doomed.size, 'pico')}` : 'Remover';
    btnConfirmRegionRemove.disabled = doomed.size === 0;
    updateApp();   // destaca a região e os picos que serão removidos
}

function closeRemoveRegionModal(skipRedraw = false) {
    if (removeModal.style.display === 'none' && !regionRemoval) return;
    removeModal.style.display = 'none';
    regionRemoval = null;
    if (!skipRedraw) updateApp();
}

[regionT0, regionT1].forEach(inp => inp.addEventListener('input', refreshRegionModal));
document.getElementById('region-select-all').addEventListener('click', () => { regionUnchecked.clear(); refreshRegionModal(); });
document.getElementById('region-select-none').addEventListener('click', () => {
    appliedPeaks.forEach(p => regionUnchecked.add(p.sample_id));
    refreshRegionModal();
});
document.getElementById('remove-region-close-btn').addEventListener('click', () => closeRemoveRegionModal());
document.getElementById('btn-cancel-region-remove').addEventListener('click', () => closeRemoveRegionModal());

btnConfirmRegionRemove.addEventListener('click', () => {
    if (!regionRemoval || regionRemoval.doomed.size === 0) return;
    const { t0, t1 } = regionRemoval;
    const set = new Set(regionRemoval.doomed);
    closeRemoveRegionModal(true);   // removePeaks já redesenha
    removePeaks(set, `Remover ${plural(set.size, 'pico')} da região ${t0.toFixed(3)}–${t1.toFixed(3)} min`);
});

// ── 2) Pela tabela: marcados ou, sem marcação, todos da guia ativa ──
// Base do "% Área". No GC-MS, na guia "Todas", todos os picos da tela são da MESMA corrida (o arquivo exibido),
// então a base é a soma de todos eles. Nos demais casos, a soma dos picos da mesma amostra/curva.
function getAreaPctBase(pk) {
    const pool = (appMode === 'gcms' && activePanelTab === 'all')
        ? domainPeaks()
        : appliedPeaks.filter(p => p.sample_id === pk.sample_id);
    return pool.reduce((acc, p) => acc + p.area, 0);
}

// No GC-FID a guia "Todas" mistura amostras (a % não faria sentido); no GC-MS é a mesma corrida
function showAreaPctColumn() {
    return !(activePanelTab === 'all' && appMode !== 'gcms');
}

function getVisiblePeaks() {
    const here = domainPeaks();
    return activePanelTab === 'all' ? here : here.filter(p => p.sample_id === activePanelTab);
}

function updateTableSelectionUI() {
    const visible = getVisiblePeaks();
    const selected = visible.filter(p => selectedPeaks.has(p));
    peaksSelectAll.checked = visible.length > 0 && selected.length === visible.length;
    peaksSelectAll.indeterminate = selected.length > 0 && selected.length < visible.length;
    btnRemovePeaks.textContent = selected.length
        ? `🗑 Remover selecionados (${selected.length})`
        : `🗑 Remover todos (${visible.length})`;
    btnRemovePeaks.disabled = visible.length === 0;
}

peaksSelectAll.addEventListener('change', () => {
    const visible = getVisiblePeaks();
    visible.forEach(p => { if (peaksSelectAll.checked) selectedPeaks.add(p); else selectedPeaks.delete(p); });
    renderBottomPanel();
});

btnRemovePeaks.addEventListener('click', () => {
    const visible = getVisiblePeaks();
    const selected = visible.filter(p => selectedPeaks.has(p));
    const target = selected.length ? selected : visible;
    if (target.length === 0) return;

    const sample = activePanelTab !== 'all' ? samples.find(s => s.id === activePanelTab) : null;
    const where = activePanelTab !== 'all' ? ` de ${sample ? sample.name : activePanelTab}` : '';
    removePeaks(new Set(target), selected.length
        ? `Remover ${plural(target.length, 'selecionado')}${where}`
        : `Remover todos os picos${where} (${target.length})`);
});

// ══════════════════════════════════════════════════════════════
// ── ALINHAR E NORMALIZAR CORRIDAS (botão ⇄) ──
// Alinhar: desloca o eixo do tempo de cada corrida para que um pico âncora caia no mesmo tR da referência.
//   Só o tempo muda (as áreas ficam iguais); picos já integrados acompanham; "Reverter" volta ao original.
// Normalizar: divide as medidas pelo mesmo pico (padrão interno) de cada corrida, ou por um fator manual.
// ══════════════════════════════════════════════════════════════
const alignBtn = document.getElementById('align-btn');
const alignModal = document.getElementById('align-modal');
const alignRef = document.getElementById('align-ref');
const alignAnchor = document.getElementById('align-anchor');
const alignWin = document.getElementById('align-win');
const alignRowsEl = document.getElementById('align-rows');
const normModeSel = document.getElementById('norm-mode');
const normIsSel = document.getElementById('norm-is');
const normMetricSel = document.getElementById('norm-metric');
const normPlotChk = document.getElementById('norm-plot');
const normRowsEl = document.getElementById('norm-rows');
const normHint = document.getElementById('norm-hint');
let alignRowsState = [];   // [{ run, apex, delta, include }]

makeDraggable(alignModal, document.getElementById('align-drag-handle'));

// Corridas do modo atual: amostras de GC-FID ou arquivos de GC-MS (sinal de referência = TIC)
function getAlignRuns() {
    if (appMode === 'gcms') return msSamples.map(s => ({ id: s.id, name: s.name, time: s.time, y: s.tic, ref: s }));
    return getFidSamples().map(s => ({ id: s.id, name: s.name, time: s.time, y: s.intensity, ref: s }));
}

function runIdOfPeak(pk) { return pk.signal ? pk.signal.ms_id : pk.sample_id; }

// Todos os picos do modo atual (não só os da tela), usados para achar o padrão interno de cada corrida
function domainAllPeaks() { return appliedPeaks.filter(pk => (appMode === 'gcms') === !!pk.signal); }

// ── Normalização ──
function getRunRef(runId) {
    const r = getAlignRuns().find(x => x.id === runId);
    return r ? r.ref : null;
}

function getNormRefPeak(runId) {
    if (!normSettings.isName) return null;
    return domainAllPeaks().find(pk => runIdOfPeak(pk) === runId && pk.peak_name === normSettings.isName) || null;
}

// Multiplicador das medidas da tabela (null = esta corrida não tem o pico padrão)
function getNormFactor(runId) {
    if (normSettings.mode === 'manual') {
        const s = getRunRef(runId);
        return s && s.normManual > 0 ? s.normManual : 1;
    }
    if (normSettings.mode === 'is') {
        const p = getNormRefPeak(runId);
        const v = p ? p[normSettings.metric] : 0;
        return v > 0 ? 1 / v : null;
    }
    return 1;
}

function getNormalizedValue(pk) {
    const f = getNormFactor(runIdOfPeak(pk));
    return f === null ? null : pk[normSettings.metric] * f;
}

function formatNormValue(pk) {
    const v = getNormalizedValue(pk);
    return v === null ? '—' : v.toLocaleString(undefined, { maximumSignificantDigits: 6 });
}

// Fator de ESCALA do cromatograma: com padrão interno, iguala a altura do padrão à da primeira corrida que o tem
function getDisplayFactor(runId) {
    if (!normSettings.plot) return 1;
    if (normSettings.mode === 'manual') {
        const s = getRunRef(runId);
        return s && s.normManual > 0 ? s.normManual : 1;
    }
    if (normSettings.mode === 'is') {
        const runs = getAlignRuns();
        const hs = runs.map(r => { const p = getNormRefPeak(r.id); return p && p.height > 0 ? p.height : null; });
        const ref = hs.find(h => h !== null);
        const mine = hs[runs.findIndex(r => r.id === runId)];
        return (ref && mine) ? ref / mine : 1;
    }
    return 1;
}

function applyDisplayFactors() {
    samples.forEach(s => { s.dispFactor = getDisplayFactor(s.ms_id !== undefined ? s.ms_id : s.id); });
}

function refreshAfterNorm() {
    updateApp();
    renderBottomPanel();
}

function renderNormPane() {
    normModeSel.value = normSettings.mode;
    normMetricSel.value = normSettings.metric;
    normPlotChk.checked = normSettings.plot;
    document.getElementById('norm-is-row').style.display = normSettings.mode === 'is' ? 'flex' : 'none';
    document.getElementById('norm-metric-row').style.display = normSettings.mode === 'none' ? 'none' : 'flex';
    document.getElementById('norm-table').parentElement.style.display = normSettings.mode === 'none' ? 'none' : 'block';
    document.getElementById('norm-th-std').textContent = normSettings.mode === 'is'
        ? `Padrão (${normSettings.metric === 'height' ? 'altura' : 'área'})` : '—';

    // candidatos a padrão interno: nomes de pico existentes (de preferência presentes em todas as corridas)
    const runs = getAlignRuns();
    const names = [...new Set(domainAllPeaks().map(p => p.peak_name))];
    names.sort((a, b) => {
        const ca = runs.filter(r => domainAllPeaks().some(p => runIdOfPeak(p) === r.id && p.peak_name === a)).length;
        const cb = runs.filter(r => domainAllPeaks().some(p => runIdOfPeak(p) === r.id && p.peak_name === b)).length;
        return cb - ca;
    });
    if (normSettings.mode === 'is' && !names.includes(normSettings.isName)) normSettings.isName = names[0] || '';
    normIsSel.innerHTML = names.length
        ? names.map(n => `<option value="${escapeHtml(n)}">${escapeHtml(n)}</option>`).join('')
        : '<option value="">(integre um pico em cada corrida)</option>';
    normIsSel.value = normSettings.isName;

    normRowsEl.innerHTML = '';
    runs.forEach(r => {
        const tr = document.createElement('tr');
        let std = '—', fac = '';
        if (normSettings.mode === 'is') {
            const p = getNormRefPeak(r.id);
            std = p ? (+p[normSettings.metric].toPrecision(6)).toLocaleString() : '<span class="warn">sem o pico</span>';
            const f = getNormFactor(r.id);
            fac = f === null ? '—' : `× ${f.toPrecision(4)}`;
        } else if (normSettings.mode === 'manual') {
            fac = `<input type="number" step="any" min="0" value="${r.ref.normManual > 0 ? r.ref.normManual : 1}">`;
        }
        tr.innerHTML = `<td>${escapeHtml(r.name)}</td><td>${std}</td><td>${fac}</td>`;
        const inp = tr.querySelector('input');
        if (inp) inp.addEventListener('change', () => {
            const v = parseFloat(inp.value);
            r.ref.normManual = v > 0 ? v : 1;
            refreshAfterNorm();
        });
        normRowsEl.appendChild(tr);
    });

    normHint.textContent = normSettings.mode === 'is'
        ? 'Coluna "Norm." da tabela = valor do pico ÷ valor do padrão na mesma corrida (sem unidade). Use o mesmo nome de pico em todas as corridas.'
        : normSettings.mode === 'manual'
            ? 'Coluna "Norm." = valor do pico × fator da corrida (por exemplo, 1 ÷ massa injetada).'
            : 'Sem normalização: a tabela mostra os valores medidos.';
}

normModeSel.addEventListener('change', () => { normSettings.mode = normModeSel.value; renderNormPane(); refreshAfterNorm(); });
normIsSel.addEventListener('change', () => { normSettings.isName = normIsSel.value; renderNormPane(); refreshAfterNorm(); });
normMetricSel.addEventListener('change', () => { normSettings.metric = normMetricSel.value; renderNormPane(); refreshAfterNorm(); });
normPlotChk.addEventListener('change', () => { normSettings.plot = normPlotChk.checked; refreshAfterNorm(); });

// ── Alinhamento ──
function apexInWindow(run, t0, t1) {
    const time = run.time, y = run.y;
    let lo = 0, hi = time.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (time[mid] < t0) lo = mid + 1; else hi = mid; }
    let best = -1;
    for (let i = lo; i < time.length && time[i] <= t1; i++) if (best < 0 || y[i] > y[best]) best = i;
    if (best < 0) return null;
    return { i: best, t: time[best], edge: best === lo || (best + 1 < time.length && time[best + 1] > t1) };
}

function populateAlignRefs() {
    const runs = getAlignRuns();
    const keep = alignRef.value;
    alignRef.innerHTML = runs.map(r => `<option value="${escapeHtml(r.id)}">${escapeHtml(r.name)}</option>`).join('');
    if (runs.some(r => r.id === keep)) alignRef.value = keep;
}

function renderAlignRows() {
    alignRowsEl.innerHTML = '';
    if (alignRowsState.length === 0) {
        alignRowsEl.innerHTML = '<tr><td colspan="5" class="muted">Informe o tR do pico âncora (ou clique em "Gráfico").</td></tr>';
        return;
    }
    alignRowsState.forEach(st => {
        const tr = document.createElement('tr');
        const shiftNow = st.run.ref.shift || 0;
        tr.innerHTML = `
            <td><input type="checkbox" ${st.include ? 'checked' : ''} ${st.apex ? '' : 'disabled'}></td>
            <td>${escapeHtml(st.run.name)}</td>
            <td>${st.apex ? st.apex.t.toFixed(4) + (st.apex.edge ? ' <span class="warn" title="O máximo está na borda da janela: aumente a janela ou ajuste a âncora">⚠</span>' : '') : '<span class="muted">não achou</span>'}</td>
            <td><input type="number" step="0.001" value="${st.delta}" ${st.apex ? '' : 'disabled'}></td>
            <td class="muted">${shiftNow ? (shiftNow > 0 ? '+' : '') + shiftNow.toFixed(4) : '—'}</td>`;
        const [cb, num] = tr.querySelectorAll('input');
        cb.addEventListener('change', () => { st.include = cb.checked; });
        num.addEventListener('change', () => { const v = parseFloat(num.value); st.delta = isNaN(v) ? 0 : v; });
        alignRowsEl.appendChild(tr);
    });
}

// Acha o ápice do pico âncora em cada corrida e calcula o deslocamento até o tR da referência
function detectAlignment() {
    const runs = getAlignRuns();
    const a = parseFloat(alignAnchor.value);
    const w = Math.max(0.005, parseFloat(alignWin.value) || 0.1);
    if (isNaN(a) || runs.length === 0) { alignRowsState = []; renderAlignRows(); return; }

    const found = runs.map(r => apexInWindow(r, a - w, a + w));
    const refIdx = runs.findIndex(r => r.id === alignRef.value);
    const target = (refIdx >= 0 && found[refIdx]) ? found[refIdx].t : a;
    alignRowsState = runs.map((r, i) => ({
        run: r, apex: found[i],
        delta: found[i] ? +(target - found[i].t).toFixed(5) : 0,
        include: !!found[i] && i !== refIdx
    }));
    renderAlignRows();
}

function shiftPeak(pk, d) {
    return {
        ...pk,
        rt_min: pk.rt_min + d,
        apex: { x: pk.apex.x + d, y: pk.apex.y },
        baseline_line: { x: pk.baseline_line.x.map(v => v + d), y: pk.baseline_line.y },
        polygon: { x: pk.polygon.x.map(v => v + d), y: pk.polygon.y }
    };
}

// Desloca UMA corrida em d minutos: tempo, picos já integrados e cache do servidor
async function shiftRun(run, d) {
    const s = run.ref, t = s.time;
    for (let i = 0; i < t.length; i++) t[i] += d;
    s.shift = Math.round(((s.shift || 0) + d) * 1e9) / 1e9;

    const isMs = appMode === 'gcms';
    for (let i = 0; i < appliedPeaks.length; i++) {
        const pk = appliedPeaks[i];
        if (isMs ? (pk.signal && pk.signal.ms_id === s.id) : (!pk.signal && pk.sample_id === s.id)) replacePeakAt(i, shiftPeak(pk, d));
    }
    try {
        if (isMs) {
            await fetch('/api/gcms/shifts', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ shifts: { [s.id]: s.shift } }) });
        } else {
            await fetch('/api/samples/shift', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ sample_id: s.id, delta: d }) });
        }
    } catch (err) {
        showToast('Aviso: o servidor não confirmou o deslocamento.');
    }
}

function afterShift() {
    resetPeakState();   // snapshots do desfazer ficariam com o tempo antigo
    updateApp();
    renderBottomPanel();
    renderAlignRows();
}

async function applyAlignment() {
    const todo = alignRowsState.filter(r => r.include && Math.abs(r.delta) > 1e-9);
    if (todo.length === 0) { showToast('Nada a deslocar: defina a âncora e marque as corridas.'); return; }
    for (const r of todo) await shiftRun(r.run, r.delta);
    showToast(`${todo.length} corrida${todo.length > 1 ? 's' : ''} alinhada${todo.length > 1 ? 's' : ''}`);
    detectAlignment();
    afterShift();
}

async function revertAlignment() {
    const runs = getAlignRuns().filter(r => Math.abs(r.ref.shift || 0) > 1e-9);
    if (runs.length === 0) { showToast('Nenhuma corrida está deslocada.'); return; }
    for (const r of runs) await shiftRun(r, -(r.ref.shift || 0));
    showToast('Alinhamento revertido: tempo original');
    detectAlignment();
    afterShift();
}

document.getElementById('btn-align-apply').addEventListener('click', applyAlignment);
document.getElementById('btn-align-revert').addEventListener('click', revertAlignment);
[alignAnchor, alignWin].forEach(inp => inp.addEventListener('change', detectAlignment));
alignRef.addEventListener('change', detectAlignment);
document.getElementById('align-pick').addEventListener('click', () => togglePicking('align-anchor'));

// ── Janela ──
function closeAlignModal() { alignModal.style.display = 'none'; }

function showAlignTab(tab) {
    alignModal.querySelectorAll('.ms-mode-tab').forEach(b => b.classList.toggle('active', b.dataset.aligntab === tab));
    document.getElementById('align-pane-align').style.display = tab === 'align' ? 'flex' : 'none';
    document.getElementById('align-pane-norm').style.display = tab === 'norm' ? 'flex' : 'none';
    if (tab === 'norm') renderNormPane();
}

alignModal.querySelectorAll('.ms-mode-tab').forEach(b => b.addEventListener('click', () => showAlignTab(b.dataset.aligntab)));

alignBtn.addEventListener('click', () => {
    if (alignModal.style.display !== 'none') { closeAlignModal(); return; }
    populateAlignRefs();
    alignRowsState = [];
    renderAlignRows();
    detectAlignment();
    renderNormPane();
    alignModal.style.display = 'flex';
});
document.getElementById('align-close-btn').addEventListener('click', closeAlignModal);

// ══════════════════════════════════════════════════════════════
// ── COMPOSTOS: lista reutilizável (método) + matriz compostos × corridas (botão 🧪) ──
// A lista vale para o modo atual (GC-FID ou GC-MS). "Aplicar" identifica/integra cada composto em cada corrida:
// reaproveita o pico já integrado na janela ou o detecta no servidor. No GC-MS cada composto pode ter o seu m/z (EIC).
// ══════════════════════════════════════════════════════════════
const compoundsBtn = document.getElementById('compounds-btn');
const compoundsModal = document.getElementById('compounds-modal');
const cmpRowsEl = document.getElementById('cmp-rows');
const cmpSnr = document.getElementById('cmp-snr');
const cmpMetric = document.getElementById('cmp-metric');
const BASELINE_TYPES = ['Valley-to-Valley', 'Drop Line', 'Tangent Skim'];

makeDraggable(compoundsModal, document.getElementById('compounds-drag-handle'));

function cmpDomain() { return appMode === 'gcms' ? 'ms' : 'fid'; }
function domainCompounds() { return compounds.filter(c => c.domain === cmpDomain()); }
function newCompoundId() { return 'c' + Date.now().toString(36) + (cmpSeq++); }

function addCompound(p = {}) {
    const c = {
        id: newCompoundId(),
        name: p.name || `Composto ${domainCompounds().length + 1}`,
        rt: p.rt !== undefined ? p.rt : null,
        win: p.win || 0.05,
        baseline: p.baseline || 'Valley-to-Valley',
        ion: p.ion !== undefined ? p.ion : null,
        tol: p.tol || msEic.tol || 0.5,
        is: false,
        domain: cmpDomain()
    };
    compounds.push(c);
    return c;
}

function renderCompoundRows() {
    const isMs = cmpDomain() === 'ms';
    document.getElementById('cmp-th-ion').style.display = isMs ? '' : 'none';
    cmpRowsEl.innerHTML = '';
    const list = domainCompounds();
    if (list.length === 0) {
        cmpRowsEl.innerHTML = `<tr><td colspan="7" class="muted">Nenhum composto. Use "+ Composto" ou "Criar dos picos".</td></tr>`;
        return;
    }
    list.forEach(c => {
        const tr = document.createElement('tr');
        tr.innerHTML = `
            <td><input type="text" class="cmp-name" value="${escapeHtml(c.name)}"></td>
            <td><input type="number" step="0.001" class="cmp-rt" value="${c.rt === null ? '' : c.rt}" placeholder="min"></td>
            <td><input type="number" step="0.01" min="0.005" class="cmp-win" value="${c.win}"></td>
            <td><select class="cmp-bl">${BASELINE_TYPES.map(b => `<option ${b === c.baseline ? 'selected' : ''}>${b}</option>`).join('')}</select></td>
            ${isMs ? `<td><input type="number" step="0.1" min="0" class="cmp-ion" value="${c.ion === null ? '' : c.ion}" placeholder="TIC"></td>` : ''}
            <td><input type="checkbox" class="cmp-is" ${c.is ? 'checked' : ''} title="Padrão interno"></td>
            <td><button class="btn-close" title="Remover composto">✕</button></td>`;
        tr.querySelector('.cmp-name').addEventListener('change', e => { c.name = e.target.value.trim() || c.name; });
        tr.querySelector('.cmp-rt').addEventListener('change', e => { const v = parseFloat(e.target.value); c.rt = isNaN(v) ? null : v; });
        tr.querySelector('.cmp-win').addEventListener('change', e => { const v = parseFloat(e.target.value); c.win = v > 0 ? v : c.win; e.target.value = c.win; });
        tr.querySelector('.cmp-bl').addEventListener('change', e => { c.baseline = e.target.value; });
        const ionIn = tr.querySelector('.cmp-ion');
        if (ionIn) ionIn.addEventListener('change', e => { const v = parseFloat(e.target.value); c.ion = v > 0 ? v : null; });
        tr.querySelector('.cmp-is').addEventListener('change', e => {
            domainCompounds().forEach(x => { x.is = false; });   // só um padrão interno por modo
            c.is = e.target.checked;
            renderCompoundRows();
        });
        tr.querySelector('.btn-close').addEventListener('click', () => {
            compounds = compounds.filter(x => x !== c);
            renderCompoundRows();
            renderCompoundMatrix();
        });
        cmpRowsEl.appendChild(tr);
    });
}

// Um composto para cada nome de pico já integrado (tR = mediana entre as corridas)
function createCompoundsFromPeaks() {
    const groups = new Map();
    domainAllPeaks().forEach(pk => {
        const key = `${pk.peak_name}|${pk.signal && pk.signal.ion !== null ? pk.signal.ion : ''}`;
        if (!groups.has(key)) groups.set(key, []);
        groups.get(key).push(pk);
    });
    let added = 0;
    groups.forEach((pks, key) => {
        const first = pks[0];
        const ion = first.signal ? first.signal.ion : null;
        if (domainCompounds().some(c => c.name === first.peak_name && c.ion === ion)) return;
        const rts = pks.map(p => p.rt_min).sort((a, b) => a - b);
        const mid = rts[Math.floor(rts.length / 2)];
        addCompound({ name: first.peak_name, rt: +mid.toFixed(4), baseline: first.baseline_type, ion, tol: first.signal ? first.signal.tol : undefined });
        added++;
    });
    renderCompoundRows();
    showToast(added ? `${added} composto${added > 1 ? 's' : ''} criado${added > 1 ? 's' : ''} dos picos integrados` : 'Nenhum nome novo nos picos integrados.');
}

function compoundSignalId(c, run) {
    return cmpDomain() === 'ms' ? msSignalId(run.id, c.ion, c.tol) : run.id;
}

async function applyCompounds() {
    const list = domainCompounds().filter(c => c.name && c.rt !== null && isFinite(c.rt) && c.win > 0);
    const runs = getAlignRuns();
    if (list.length === 0) { showToast('Defina ao menos um composto com nome e tR.'); return; }
    if (runs.length === 0) { showToast('Não há corridas carregadas.'); return; }
    const snr = Math.max(0, parseFloat(cmpSnr.value) || 0);

    const pairs = [];
    list.forEach(c => runs.forEach(r => pairs.push({ c, r, sid: compoundSignalId(c, r) })));

    // 1) picos já integrados na janela são reaproveitados; 2) os demais vão para a detecção no servidor
    const reuse = [], items = [];
    pairs.forEach(pr => {
        const cands = appliedPeaks
            .filter(pk => pk.sample_id === pr.sid && Math.abs(pk.rt_min - pr.c.rt) <= pr.c.win && (!pk.compound_id || pk.compound_id === pr.c.id))
            .sort((a, b) => Math.abs(a.rt_min - pr.c.rt) - Math.abs(b.rt_min - pr.c.rt));
        if (cands.length) reuse.push({ pk: cands[0], pr });
        else items.push({ signal_id: pr.sid, compound_id: pr.c.id, rt: pr.c.rt, win: pr.c.win, baseline_type: pr.c.baseline });
    });

    let results = [];
    if (items.length) {
        try {
            const res = await fetch('/api/compounds/detect', {
                method: 'POST', headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ items, snr_min: snr })
            });
            if (!res.ok) throw new Error('O servidor recusou o pedido.');
            results = (await res.json()).results;
        } catch (err) {
            showToast('Falha ao detectar picos: ' + err.message);
            return;
        }
    }

    pushHistory(`Aplicar método de compostos (${list.length} × ${runs.length})`);

    // resultados antigos destes compostos (que não foram reaproveitados) saem
    const keep = new Set(reuse.map(x => x.pk));
    const affected = new Set(pairs.map(p => `${p.c.id}|${p.sid}`));
    appliedPeaks = appliedPeaks.filter(pk => !(pk.compound_id && affected.has(`${pk.compound_id}|${pk.sample_id}`) && !keep.has(pk)));

    reuse.forEach(({ pk, pr }) => {
        const idx = appliedPeaks.indexOf(pk);
        if (idx > -1) replacePeakAt(idx, { ...pk, peak_name: pr.c.name, compound_id: pr.c.id });
    });

    const byId = Object.fromEntries(list.map(c => [c.id, c]));
    let detected = 0;
    results.forEach(r => {
        if (!r.found) return;
        const c = byId[r.compound_id];
        appliedPeaks.push({ ...r.peak, peak_name: c.name, compound_id: c.id });
        detected++;
    });
    const total = pairs.length;
    const nd = total - reuse.length - detected;

    selectedPeaks.forEach(p => { if (!appliedPeaks.includes(p)) selectedPeaks.delete(p); });
    updateApp();
    renderBottomPanel();
    renderCompoundMatrix();
    showCmpTab('matrix');
    showToast(`Método aplicado: ${reuse.length} reaproveitado${reuse.length === 1 ? '' : 's'}, ${detected} detectado${detected === 1 ? '' : 's'}, ${nd} não detectado${nd === 1 ? '' : 's'}`);
}

// ── Matriz ──
function getCompoundPeak(c, runId) {
    return appliedPeaks.find(pk => pk.compound_id === c.id && runIdOfPeak(pk) === runId && (!!pk.signal) === (c.domain === 'ms')) || null;
}

function cmpCellValue(c, run, metric, list) {
    const pk = getCompoundPeak(c, run.id);
    if (!pk) return null;
    if (metric === 'area') return pk.area;
    if (metric === 'height') return pk.height;
    if (metric === 'rt') return pk.rt_min;
    if (metric === 'pct') {
        const tot = list.reduce((acc, x) => { const p = getCompoundPeak(x, run.id); return acc + (p ? p.area : 0); }, 0);
        return tot > 0 ? (pk.area / tot) * 100 : null;
    }
    if (metric === 'ratio') {
        const is = list.find(x => x.is);
        const pis = is ? getCompoundPeak(is, run.id) : null;
        return pis && pis.area > 0 ? pk.area / pis.area : null;
    }
    return null;
}

function formatCmpValue(v, metric) {
    if (v === null) return 'n.d.';
    if (metric === 'rt') return v.toFixed(4);
    if (metric === 'pct') return v.toFixed(2) + ' %';
    if (metric === 'ratio') return v.toLocaleString(undefined, { maximumSignificantDigits: 5 });
    return v.toLocaleString(undefined, { minimumFractionDigits: 2, maximumFractionDigits: 2 });
}

function renderCompoundMatrix() {
    const list = domainCompounds().filter(c => c.name);
    const runs = getAlignRuns();
    const metric = cmpMetric.value;
    const head = document.getElementById('cmp-matrix-head');
    const body = document.getElementById('cmp-matrix-body');
    head.innerHTML = `<tr><th>Composto</th><th>tR esp.</th>${runs.map(r => `<th>${escapeHtml(r.name)}</th>`).join('')}</tr>`;
    body.innerHTML = '';
    if (list.length === 0 || runs.length === 0) {
        body.innerHTML = `<tr><td colspan="${2 + runs.length}" class="muted">Defina compostos na aba Lista e aplique às corridas.</td></tr>`;
        return;
    }
    list.forEach(c => {
        const tr = document.createElement('tr');
        let cells = `<td>${escapeHtml(c.name)}${c.is ? ' <span class="muted" title="Padrão interno">(PI)</span>' : ''}${c.ion !== null ? ` <span class="muted">m/z ${+c.ion}</span>` : ''}</td><td class="num">${c.rt === null ? '—' : c.rt.toFixed(3)}</td>`;
        runs.forEach(r => {
            const v = cmpCellValue(c, r, metric, list);
            cells += `<td class="${v === null ? 'nd' : 'num'}">${formatCmpValue(v, metric)}</td>`;
        });
        tr.innerHTML = cells;
        body.appendChild(tr);
    });
}

function exportCompoundMatrix() {
    const list = domainCompounds().filter(c => c.name);
    const runs = getAlignRuns();
    if (list.length === 0 || runs.length === 0) { alert('Não há compostos/corridas para exportar.'); return; }
    const metric = cmpMetric.value;
    const label = { area: 'Area (u.a.*s)', height: 'Altura', rt: 'tR (min)', pct: '% Area', ratio: 'Area / padrao interno' }[metric];
    let csv = '\uFEFF' + `# ${label}\n` + ['Composto', 'tR esperado (min)', ...runs.map(r => r.name)].join(';') + '\n';
    list.forEach(c => {
        const row = [`"${c.name}${c.ion !== null ? ' (m/z ' + c.ion + ')' : ''}"`, c.rt === null ? '' : c.rt];
        runs.forEach(r => { const v = cmpCellValue(c, r, metric, list); row.push(v === null ? 'n.d.' : String(v)); });
        csv += row.join(';') + '\n';
    });
    downloadText(csv, `matriz_compostos_${metric}.csv`, 'text/csv;charset=utf-8;');
}

// ── Janela ──
function showCmpTab(tab) {
    compoundsModal.querySelectorAll('.ms-mode-tab').forEach(b => b.classList.toggle('active', b.dataset.cmptab === tab));
    document.getElementById('cmp-pane-list').style.display = tab === 'list' ? 'flex' : 'none';
    document.getElementById('cmp-pane-matrix').style.display = tab === 'matrix' ? 'flex' : 'none';
    if (tab === 'matrix') renderCompoundMatrix();
}

function closeCompoundsModal() { compoundsModal.style.display = 'none'; }

compoundsBtn.addEventListener('click', () => {
    if (compoundsModal.style.display !== 'none') { closeCompoundsModal(); return; }
    renderCompoundRows();
    renderCompoundMatrix();
    compoundsModal.style.display = 'flex';
});
compoundsModal.querySelectorAll('.ms-mode-tab').forEach(b => b.addEventListener('click', () => showCmpTab(b.dataset.cmptab)));
document.getElementById('compounds-close-btn').addEventListener('click', closeCompoundsModal);
document.getElementById('cmp-add').addEventListener('click', () => { addCompound(); renderCompoundRows(); });
document.getElementById('cmp-from-peaks').addEventListener('click', createCompoundsFromPeaks);
document.getElementById('cmp-apply').addEventListener('click', applyCompounds);
document.getElementById('cmp-export-csv').addEventListener('click', exportCompoundMatrix);
cmpMetric.addEventListener('change', renderCompoundMatrix);
document.getElementById('cmp-clear-peaks').addEventListener('click', () => {
    const ids = new Set(domainCompounds().map(c => c.id));
    const doomed = new Set(appliedPeaks.filter(pk => pk.compound_id && ids.has(pk.compound_id) && (appMode === 'gcms') === !!pk.signal));
    if (doomed.size === 0) { showToast('Não há resultados de compostos para limpar.'); return; }
    removePeaks(doomed, `Limpar resultados do método (${doomed.size})`);
    renderCompoundMatrix();
});

// ══════════════════════════════════════════════════════════════
// ── GRUPOS DE AMOSTRAS (só organização) ──
// Arraste uma amostra sobre outra para criar um grupo; sobre um grupo para entrar nele; para o espaço vazio da lista
// para sair. No GC-FID, clicar no grupo mostra/esmaece todas as amostras dele no gráfico. Duplo clique no nome renomeia.
// ══════════════════════════════════════════════════════════════
let dragSample = null;   // { id, domain } enquanto uma amostra está sendo arrastada
let groupSeq = 1;

function getListItems(domain) {
    return domain === 'ms' ? msSamples : getFidSamples();
}

// Remove grupos sem membros e referências a grupos que não existem mais
function pruneGroups() {
    sampleGroups = sampleGroups.filter(g => getListItems(g.domain).some(s => s.group === g.id));
    ['fid', 'ms'].forEach(domain => {
        getListItems(domain).forEach(s => {
            if (s.group && !sampleGroups.some(g => g.id === s.group && g.domain === domain)) delete s.group;
        });
    });
}

function newGroupFor(domain) {
    let n = sampleGroups.filter(g => g.domain === domain).length + 1;
    while (sampleGroups.some(g => g.domain === domain && g.name === `Grupo ${n}`)) n++;
    const g = { id: 'g' + Date.now().toString(36) + (groupSeq++), name: `Grupo ${n}`, collapsed: false, domain };
    sampleGroups.push(g);
    return g;
}

// Remove os realces de arrastar. Necessário porque a lista é redesenhada durante o drop e o navegador
// não dispara "dragend" em elementos que saíram da tela.
function clearDragMarks() {
    dragSample = null;
    document.querySelectorAll('.drop-target, .dragging, .drop-ungroup').forEach(el => el.classList.remove('drop-target', 'dragging', 'drop-ungroup'));
}

function afterGroupChange() {
    clearDragMarks();
    pruneGroups();
    renderSampleList();
}

function findDragged(drag) {
    return drag ? getListItems(drag.domain).find(s => s.id === drag.id) : null;
}

function dropOnSample(drag, target, domain) {
    const a = findDragged(drag);
    if (!a || drag.domain !== domain || a === target) return;
    if (target.group) {
        if (a.group === target.group) return;
        a.group = target.group;
    } else {
        const g = newGroupFor(domain);   // arrastar sobre uma amostra solta cria um grupo novo com as duas
        target.group = g.id;
        a.group = g.id;
    }
    afterGroupChange();
}

function dropOnGroup(drag, g) {
    const a = findDragged(drag);
    if (!a || drag.domain !== g.domain || a.group === g.id) return;
    a.group = g.id;
    afterGroupChange();
}

function ungroupSample(drag) {
    const a = findDragged(drag);
    if (!a || !a.group) return;
    delete a.group;
    afterGroupChange();
}

function attachSampleDnD(li, s, domain) {
    li.draggable = true;
    li.addEventListener('dragstart', (e) => {
        dragSample = { id: s.id, domain };
        e.dataTransfer.effectAllowed = 'move';
        e.dataTransfer.setData('text/plain', s.id);
        li.classList.add('dragging');
    });
    li.addEventListener('dragend', () => {
        dragSample = null;
        document.querySelectorAll('#sample-list .dragging, #sample-list .drop-target').forEach(el => el.classList.remove('dragging', 'drop-target'));
        sampleList.parentElement.classList.remove('drop-ungroup');
    });
    li.addEventListener('dragover', (e) => {
        if (!dragSample || dragSample.domain !== domain || dragSample.id === s.id) return;
        e.preventDefault();
        li.classList.add('drop-target');
    });
    li.addEventListener('dragleave', () => li.classList.remove('drop-target'));
    li.addEventListener('drop', (e) => {
        e.preventDefault();
        e.stopPropagation();
        li.classList.remove('drop-target');
        dropOnSample(dragSample, s, domain);
        clearDragMarks();
    });
    return li;
}

function makeGroupHeader(g, members, domain) {
    const li = document.createElement('li');
    const allHidden = domain === 'fid' && members.every(m => !m.visible);
    li.className = 'sample-item sample-group' + (allHidden ? ' is-hidden' : '');
    li.title = domain === 'fid'
        ? 'Clique: mostrar / esmaecer o grupo no gráfico · Duplo clique no nome: renomear'
        : 'Clique: recolher / expandir · Duplo clique no nome: renomear';
    li.innerHTML = `
        <div class="sample-left">
            <span class="group-caret" title="Recolher / expandir">${g.collapsed ? '▸' : '▾'}</span>
            <span class="sample-name group-name">${escapeHtml(g.name)}</span>
            <span class="group-count">${members.length}</span>
        </div>
        <button class="btn-close" title="Desfazer o grupo (as amostras continuam)">✕</button>`;

    // A ação do clique simples é adiada: se virar duplo clique (renomear), ela é cancelada.
    // (Sem isso, o 1º clique redesenha a lista e o duplo clique nunca chega ao nome.)
    let clickTimer = null;
    li.addEventListener('click', (e) => {
        if (e.target.closest('.rename-input')) return;
        if (e.target.closest('.group-caret')) {
            g.collapsed = !g.collapsed;
            renderSampleList();
            return;
        }
        if (e.detail > 1) return;   // 2º clique de um duplo clique
        clearTimeout(clickTimer);
        clickTimer = setTimeout(() => {
            if (domain === 'ms') {
                g.collapsed = !g.collapsed;
                renderSampleList();
                return;
            }
            const anyVisible = members.some(m => m.visible);   // esconde se algum está visível, senão mostra
            members.forEach(m => { m.visible = !anyVisible; });
            updateApp();
            renderBottomPanel();
        }, 230);
    });

    li.querySelector('.group-name').addEventListener('dblclick', (e) => {
        e.stopPropagation();
        clearTimeout(clickTimer);
        const nameSpan = e.currentTarget;
        const current = g.name;
        const input = document.createElement('input');
        input.type = 'text';
        input.className = 'rename-input';
        input.value = current;
        li.draggable = false;
        nameSpan.replaceWith(input);
        input.focus();
        input.select();
        let cancelled = false;
        const finish = () => { if (!cancelled) g.name = input.value.trim() || current; renderSampleList(); };
        input.addEventListener('keydown', (ke) => { if (ke.key === 'Enter') finish(); if (ke.key === 'Escape') { cancelled = true; renderSampleList(); } });
        input.addEventListener('blur', finish);
    });

    li.querySelector('.btn-close').addEventListener('click', (e) => {
        e.stopPropagation();
        members.forEach(m => { delete m.group; });
        sampleGroups = sampleGroups.filter(x => x !== g);
        renderSampleList();
    });

    li.addEventListener('dragover', (e) => {
        if (!dragSample || dragSample.domain !== domain) return;
        e.preventDefault();
        li.classList.add('drop-target');
    });
    li.addEventListener('dragleave', () => li.classList.remove('drop-target'));
    li.addEventListener('drop', (e) => {
        e.preventDefault();
        e.stopPropagation();
        li.classList.remove('drop-target');
        dropOnGroup(dragSample, g);
        clearDragMarks();
    });
    return li;
}

// Desenha a lista: amostras soltas e grupos (cabeçalho + membros), na ordem da primeira aparição
function renderGroupedList(items, domain, makeLi) {
    pruneGroups();
    sampleList.innerHTML = '';
    const done = new Set();
    items.forEach(s => {
        const g = s.group ? sampleGroups.find(x => x.id === s.group && x.domain === domain) : null;
        if (!g) {
            sampleList.appendChild(attachSampleDnD(makeLi(s, false), s, domain));
            return;
        }
        if (done.has(g.id)) return;
        done.add(g.id);
        const members = items.filter(m => m.group === g.id);
        sampleList.appendChild(makeGroupHeader(g, members, domain));
        if (!g.collapsed) members.forEach(m => sampleList.appendChild(attachSampleDnD(makeLi(m, true), m, domain)));
    });
}

// Soltar no espaço vazio da lista tira a amostra do grupo
(() => {
    const area = sampleList.parentElement;
    area.addEventListener('dragover', (e) => {
        if (!dragSample || e.target.closest('li')) return;
        e.preventDefault();
        area.classList.add('drop-ungroup');
    });
    area.addEventListener('dragleave', (e) => { if (!area.contains(e.relatedTarget)) area.classList.remove('drop-ungroup'); });
    area.addEventListener('drop', (e) => {
        if (!dragSample || e.target.closest('li')) return;
        e.preventDefault();
        area.classList.remove('drop-ungroup');
        ungroupSample(dragSample);
        clearDragMarks();
    });
})();

// ══════════════════════════════════════════════════════════════
// ── MODO GC-MS ──
// Vários .CDF por projeto, mas só um é exibido por vez. O gráfico mostra o TIC;
// a gaveta inferior mostra o espectro de massas (fragmentação) de um scan.
// Duplo clique no cromatograma (ou ← →) muda o scan exibido.
// ══════════════════════════════════════════════════════════════
const modeBtns = document.querySelectorAll('.mode-btn');
const msPanel = document.getElementById('ms-panel');
const msPanelResizer = document.getElementById('ms-panel-resizer');
const msTitle = document.getElementById('ms-title');
const msSpectrumDiv = document.getElementById('ms-spectrum');
const msEmptyMsg = document.getElementById('ms-empty-msg');
const msRelBtn = document.getElementById('ms-rel-btn');
const msExportBtn = document.getElementById('ms-export-btn');
const msMinimizeBtn = document.getElementById('ms-minimize-btn');
const msExportMspBtn = document.getElementById('ms-export-msp-btn');
const msAvgBtn = document.getElementById('ms-avg-btn');
const msModeTabs = document.getElementById('ms-mode-tabs');
const msNav = document.getElementById('ms-nav');
const msAvgControls = document.getElementById('ms-avg-controls');
const msBinSelect = document.getElementById('ms-bin-select');
const msBgBtn = document.getElementById('ms-bg-btn');
const msBgToggle = document.getElementById('ms-bg-toggle');

const msSpectrumCache = new Map();
let msSpectrumToken = 0;

function escapeHtml(str) {
    return String(str).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

function getActiveMs() {
    return msSamples.find(m => m.id === activeMsId) || null;
}

function clampScan(s, idx) {
    return Math.max(0, Math.min(s.time.length - 1, Math.round(idx)));
}

// Índice do scan mais próximo de um tempo (busca binária; time é crescente)
function findClosestTimeIndex(times, t) {
    let lo = 0, hi = times.length - 1;
    while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (times[mid] < t) lo = mid; else hi = mid;
    }
    return Math.abs(times[lo] - t) <= Math.abs(times[hi] - t) ? lo : hi;
}

function getMsMarkerShapes() {
    const s = getActiveMs();
    if (!s || !s.time.length) return [];
    const shapes = [];
    const band = (r, fill, border) => ({
        type: 'rect', xref: 'x', yref: 'paper',
        x0: s.time[r.i0], x1: s.time[r.i1], y0: 0, y1: 1,
        fillcolor: fill, line: { color: border, width: 1, dash: 'dot' }, layer: 'below'
    });

    // intervalo do espectro médio (azul) e região de fundo (âmbar)
    if (s.avg) {
        shapes.push(band(s.avg, 'rgba(2, 132, 199, 0.13)', 'rgba(2, 132, 199, 0.75)'));
        if (s.avg.bg) shapes.push(band(s.avg.bg, 'rgba(245, 158, 11, 0.15)', 'rgba(245, 158, 11, 0.85)'));
    }
    // posição do scan exibido (modo Scan)
    if (msViewMode() === 'scan') {
        const t = s.time[clampScan(s, s.scanIdx)];
        shapes.push({
            type: 'line', xref: 'x', yref: 'paper', x0: t, x1: t, y0: 0, y1: 1,
            line: { color: '#f59e0b', width: 1.5, dash: 'dash' }
        });
    }
    // pré-visualização do início escolhido durante a seleção
    if (pickState === 'end' && (pickPurpose === 'ms-avg' || pickPurpose === 'ms-bg') && pickedStartTime !== null) {
        shapes.push({
            type: 'line', xref: 'x', yref: 'paper', x0: pickedStartTime, x1: pickedStartTime, y0: 0, y1: 1,
            line: { color: pickPurpose === 'ms-bg' ? '#f59e0b' : '#0284c7', width: 1.5, dash: 'dot' }
        });
    }
    return shapes;
}

// relayout de autorange pode falhar dentro do Plotly quando só há a rampa (eixo y2) no gráfico
function safeRelayout(update) {
    try { Plotly.relayout(plotDiv, update); } catch (err) { /* gráfico sem dados: nada a ajustar */ }
}

function msAutorange() {
    const u = { 'xaxis.autorange': true };
    if (isYPercentMode) { u['yaxis.range'] = [0, 105]; u['yaxis.autorange'] = false; }
    else u['yaxis.autorange'] = true;
    safeRelayout(u);
}

// ── Troca de modo ──
function setMode(mode, force = false) {
    if (!force && mode === appMode) return;
    const leavingFid = appMode === 'gcfid' && mode === 'gcms';
    if (leavingFid) {
        if (pickState !== 'idle') cancelPicking();
        if (modal.style.display !== 'none') { modal.style.display = 'none'; setIntegrationModalMode(false); }
        closeRemoveRegionModal(true);
    }
    appMode = mode;
    activateDomain();   // `samples` = amostras do FID ou sinais do GC-MS
    const isMs = mode === 'gcms';

    modeBtns.forEach(b => b.classList.toggle('active', b.dataset.mode === mode));
    msAvgBtn.style.display = isMs ? '' : 'none';
    eicBtn.style.display = isMs ? '' : 'none';
    contextMenu.classList.toggle('ms-mode', isMs);

    if (alignModal.style.display !== 'none') { populateAlignRefs(); detectAlignment(); renderNormPane(); }
    if (compoundsModal.style.display !== 'none') { renderCompoundRows(); renderCompoundMatrix(); }
    updateApp();
    if (isMs) {
        showMsPanelIfNeeded();
        renderBottomPanel();
        msAutorange();
    } else {
        if (pickState !== 'idle') cancelPicking();
        closeEicModal();
        hideMsPanel();
        renderBottomPanel();
        safeRelayout({
            'xaxis.autorange': true,
            'yaxis.autorange': !isYPercentMode,
            ...(isYPercentMode ? { 'yaxis.range': [0, 105] } : {})
        });
    }
    setTimeout(() => Plotly.Plots.resize(plotDiv), 60);
}

modeBtns.forEach(b => b.addEventListener('click', () => setMode(b.dataset.mode)));

function clearMsState() {
    msSamples = [];
    activeMsId = null;
    msCurrentSpectrum = null;
    msAvgResult = null;
    msPanelMode = 'scan';
    msReference = null;
    msCmpResult = null;
    msCmpSource = 'scan';
    msEic = { ions: [], tol: 0.5, showTic: false };
    msEicData = null;
    msEicToken++;
    syncEicInputs();
    closeEicModal();
    msSpectrumCache.clear();
    msSpectrumToken++;
    hideMsPanel();
    fetch('/api/gcms/release', { method: 'POST' }).catch(() => { });
}

// ── Cromatograma (TIC) ──
function updateMsApp() {
    rebuildMsSignals();
    applyDisplayFactors();
    sampleCount.innerText = msSamples.length;
    const isDark = currentTheme === 'dark';
    const plotData = [];
    const s = getActiveMs();
    const f = (s && !isYPercentMode) ? (samples[0] ? (samples[0].dispFactor || 1) : 1) : 1;   // normalização (só em escala absoluta)

    if (s) {
        if (!s.maxTic) s.maxTic = s.tic.reduce((a, b) => (b > a ? b : a), 0) || 1;
        const i = clampScan(s, s.scanIdx);

        if (isEicActive()) {
            // EIC: uma curva por íon. A escala % usa o maior valor entre todos os íons (preserva as proporções).
            const series = msEicData.series;
            let maxE = 0;
            series.forEach(arr => { const m = arrayMax(arr); if (m > maxE) maxE = m; });
            maxE = maxE || 1;
            msEicData.signals.forEach(sig => { sig.maxVal = maxE; });   // mesma escala % das curvas

            if (msEic.showTic) {
                // TIC ao fundo, reescalado para ter a mesma altura máxima das curvas (só referência visual)
                const k = isYPercentMode ? 100 / s.maxTic : (maxE * f) / s.maxTic;
                plotData.push({
                    x: s.time, y: s.tic.map(v => v * k), name: 'TIC (escala ajustada)',
                    type: 'scattergl', mode: 'lines',
                    line: { color: isDark ? 'rgba(148,163,184,0.45)' : 'rgba(100,116,139,0.40)', width: 1 },
                    hoverinfo: 'none'
                });
            }
            series.forEach((arr, k) => {
                plotData.push({
                    x: s.time, y: isYPercentMode ? arr.map(v => (v / maxE) * 100) : (f === 1 ? arr : arr.map(v => v * f)),
                    name: `m/z ${+msEic.ions[k].toFixed(2)}`,
                    type: 'scattergl', mode: 'lines',
                    line: { color: EIC_COLORS[k % EIC_COLORS.length], width: 1.4 },
                    hoverinfo: 'none'
                });
            });
        } else {
            if (s._ticSig) s._ticSig.maxVal = s.maxTic;
            const y = isYPercentMode ? s.tic.map(v => (v / s.maxTic) * 100) : (f === 1 ? s.tic : s.tic.map(v => v * f));
            plotData.push({
                x: s.time, y: y, name: s.name,
                type: 'scattergl', mode: 'lines',
                line: { color: s.color, width: 1.2 },
                hoverinfo: 'none'
            });
            if (msViewMode() === 'scan') {
                plotData.push({
                    x: [s.time[i]], y: [y[i]],
                    type: 'scatter', mode: 'markers',
                    marker: { size: 8, color: '#f59e0b', line: { color: isDark ? '#f8fafc' : '#0f172a', width: 1.5 } },
                    hoverinfo: 'none'
                });
            }
        }
    }

    // Picos integrados e prévia da integração (mesmas rotinas do GC-FID)
    getAppliedPeakTraces().concat(getIntegrationPreviewTraces()).forEach(tr => { tr.showlegend = false; plotData.push(tr); });

    // Rampas de temperatura (mantidas do GC-FID)
    const visibleRamps = (temperatureRamps || []).filter(r => r.visible);
    if (visibleRamps.length > 0) {
        const maxPlotTime = s && s.time.length ? s.time[s.time.length - 1] : 0;
        visibleRamps.forEach(ramp => {
            const pts = calculateRampPoints(ramp, maxPlotTime);
            if (pts.x && pts.x.length > 0) {
                plotData.push({
                    x: pts.x, y: pts.y, name: ramp.name || 'Rampa de Temp',
                    type: 'scatter', mode: 'lines',
                    line: { dash: 'dot', width: 1.8, color: isDark ? '#64748b' : '#94a3b8' },
                    yaxis: 'y2', hoverinfo: 'none'
                });
            }
        });
    }

    Plotly.react(plotDiv, plotData, getPlotLayout(), getPlotConfig());
    renderMsSampleList();
}

// ── Lista de amostras (GC-MS) ──
function selectMsSample(id) {
    if (id === activeMsId) return;
    activeMsId = id;
    msCurrentSpectrum = null;
    msAvgResult = null;
    updateApp();
    showMsPanelIfNeeded();
    msAutorange();
}

function renderMsSampleList() {
    renderGroupedList(msSamples, 'ms', (s, inGroup) => {
        const li = document.createElement('li');
        li.className = `sample-item ms-sample ${s.id === activeMsId ? 'ms-active' : ''} ${inGroup ? 'in-group' : ''}`;
        li.title = 'Clique para exibir este arquivo';
        li.innerHTML = `
            <div class="sample-left">
                <div class="color-picker-wrapper" title="Mudar cor">
                    <span class="color-dot" style="background-color: ${s.color}"></span>
                    <input type="color" class="color-input-hidden" value="${s.color}">
                </div>
                <span class="sample-name" title="Duplo clique p/ renomear">${escapeHtml(s.name)}</span>
            </div>
            <button class="btn-close" title="Remover arquivo">✕</button>
        `;

        let clickTimer = null;   // clique adiado: não atrapalha o duplo clique de renomear
        li.addEventListener('click', (e) => {
            if (e.target.closest('.color-picker-wrapper') || e.target.closest('.rename-input')) return;
            if (e.detail > 1) return;
            clearTimeout(clickTimer);
            clickTimer = setTimeout(() => selectMsSample(s.id), 230);
        });

        const colorWrapper = li.querySelector('.color-picker-wrapper');
        const colorInput = li.querySelector('.color-input-hidden');
        colorWrapper.addEventListener('click', (e) => { e.stopPropagation(); colorInput.click(); });
        colorInput.addEventListener('change', (e) => { s.color = e.target.value; updateApp(); renderSpectrum(); });

        const nameSpan = li.querySelector('.sample-name');
        nameSpan.addEventListener('dblclick', (e) => {
            e.stopPropagation();
            clearTimeout(clickTimer);
            const currentName = s.name;
            const input = document.createElement('input');
            input.type = 'text';
            input.className = 'rename-input';
            input.value = currentName;
            nameSpan.replaceWith(input);
            input.focus();
            input.select();
            const finish = () => { s.name = input.value.trim() || currentName; updateApp(); updateMsTitle(); };
            input.addEventListener('keydown', (ke) => { if (ke.key === 'Enter') finish(); if (ke.key === 'Escape') updateApp(); });
            input.addEventListener('blur', finish);
        });

        li.querySelector('.btn-close').addEventListener('click', (e) => {
            e.stopPropagation();
            msSamples = msSamples.filter(item => item.id !== s.id);
            appliedPeaks = appliedPeaks.filter(p => !(p.signal && p.signal.ms_id === s.id));
            resetPeakState();
            if (activeMsId === s.id) {
                activeMsId = msSamples[0] ? msSamples[0].id : null;
                msCurrentSpectrum = null;
                msAvgResult = null;
            }
            updateApp();
            showMsPanelIfNeeded();
            msAutorange();
        });

        return li;
    });
}

// ── Upload de .CDF de GC-MS (vários por vez) ──
async function uploadMsFiles(files) {
    const formData = new FormData();
    for (let i = 0; i < files.length; i++) formData.append('files', files[i]);

    uploadBtn.disabled = true;
    try {
        const res = await fetch('/api/upload/gc-ms', { method: 'POST', body: formData });
        const data = await res.json();

        if (data.errors && data.errors.length > 0) {
            alert('Não foi possível ler:\n' + data.errors.map(e => `• ${e.filename}: ${e.error}`).join('\n'));
        }

        if (data.success && data.success.length > 0) {
            data.success.forEach(m => {
                const existing = msSamples.findIndex(x => x.id === m.sample_id);
                const entry = {
                    id: m.sample_id,
                    name: sanitizeFilename(m.filename),
                    color: existing >= 0 ? msSamples[existing].color : COLOR_PALETTE[msSamples.length % COLOR_PALETTE.length],
                    n_scans: m.n_scans,
                    mz_min: m.mz_min,
                    mz_max: m.mz_max,
                    time: m.time,
                    tic: m.tic,
                    scanIdx: Math.floor(m.time.length / 2)   // começa no meio do cromatograma
                };
                if (existing >= 0) msSamples[existing] = entry; else msSamples.push(entry);
                if (!getActiveMs()) activeMsId = entry.id;   // o primeiro arquivo importado já é exibido
            });
            msSpectrumCache.clear();
            msCurrentSpectrum = null;
            msAvgResult = null;
            msEicData = null;
            updateApp();
            showMsPanelIfNeeded();
            msAutorange();
        }
    } catch (err) {
        alert('Erro ao carregar arquivos GC-MS.');
    } finally {
        uploadBtn.disabled = false;
        fileInput.value = '';
    }
}

// ── Navegação entre scans ──
function goToScan(s, idx) {
    s.scanIdx = clampScan(s, idx);
    if (msPanelMode === 'cmp') msCmpSource = 'scan';   // comparando: o novo scan passa a ser o comparado
    else msPanelMode = 'scan';                         // senão, escolher um scan leva a gaveta ao modo Scan
    updateMsPanelUI();
    updateApp();
    loadSpectrum(s, s.scanIdx);
}

function stepScan(delta) {
    const s = getActiveMs();
    if (s) goToScan(s, s.scanIdx + delta);
}

// Duplo clique no cromatograma: leva o espectro para aquele tempo.
// O Plotly intercepta o mouse e o evento nativo 'dblclick' não chega aqui, então o duplo clique
// é detectado a partir de dois 'click' seguidos (< 450 ms e quase na mesma posição).
let msLastClick = { t: 0, x: 0, y: 0 };
plotDiv.addEventListener('click', (e) => {
    if (appMode !== 'gcms' || pickState !== 'idle') return;   // durante a seleção de intervalo, cliques não escolhem scan
    const now = performance.now();
    const isDouble = (now - msLastClick.t) < 450 && Math.hypot(e.clientX - msLastClick.x, e.clientY - msLastClick.y) < 8;
    msLastClick = isDouble ? { t: 0, x: 0, y: 0 } : { t: now, x: e.clientX, y: e.clientY };
    if (!isDouble) return;

    const s = getActiveMs();
    const fl = plotDiv._fullLayout;
    if (!s || !fl || !fl.xaxis) return;

    const rect = plotDiv.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const xa = fl.xaxis;
    if (mouseX < xa._offset || mouseX > xa._offset + xa._length) return;

    goToScan(s, findClosestTimeIndex(s.time, xa.p2c(mouseX - xa._offset)));
});

[['ms-step-back-10', -10], ['ms-step-back-1', -1], ['ms-step-fwd-1', 1], ['ms-step-fwd-10', 10]].forEach(([id, d]) => {
    document.getElementById(id).addEventListener('click', () => stepScan(d));
});

document.addEventListener('keydown', (e) => {
    if (appMode !== 'gcms') return;
    if (e.ctrlKey || e.metaKey || e.altKey) return;
    const tag = (e.target.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'textarea' || tag === 'select' || e.target.isContentEditable) return;
    if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
        e.preventDefault();
        stepScan((e.key === 'ArrowRight' ? 1 : -1) * (e.shiftKey ? 10 : 1));
    }
});

// ── Gaveta do espectro de massas ──
function hideMsPanel() {
    if (msPanel) msPanel.style.display = 'none';
    setTimeout(() => Plotly.Plots.resize(plotDiv), 60);
}

function showMsPanelIfNeeded() {
    const s = getActiveMs();
    if (appMode !== 'gcms' || !s) { hideMsPanel(); return; }
    if (msPanelMode === 'avg' && !s.avg) msPanelMode = 'scan';
    if (msPanelMode === 'cmp' && !msReference) msPanelMode = 'scan';
    if (msEic.ions.length > 0 && !isEicActive()) refreshMsEic();   // outra amostra: recalcula os mesmos íons
    msPanel.style.display = 'flex';
    updateMsPanelUI();
    renderBottomPanel();   // a tabela mostra os picos dos sinais da amostra ativa
    if (msViewMode() === 'avg') refreshMsAvg(); else loadSpectrum(s, s.scanIdx);
    setTimeout(() => Plotly.Plots.resize(plotDiv), 60);
}

function formatBig(v) {
    const a = Math.abs(v);
    if (a >= 1e9) return (v / 1e9).toFixed(2) + 'G';
    if (a >= 1e6) return (v / 1e6).toFixed(2) + 'M';
    if (a >= 1e3) return (v / 1e3).toFixed(1) + 'k';
    return v.toFixed(0);
}

function updateMsTitle() {
    const s = getActiveMs();
    if (!s) return;
    const dim = (txt) => ` <span class="ms-dim">${txt}</span>`;
    const spec = getDisplaySpectrum();
    let html;

    if (msPanelMode === 'cmp' && msReference) {
        const r = msCmpResult;
        const what = msViewMode() === 'avg' ? 'média' : `scan ${clampScan(s, s.scanIdx) + 1}`;
        html = `<strong>Comparação</strong>${dim('·')} ${escapeHtml(s.name)} ${dim(`(${what})`)} × ${escapeHtml(msReference.name)}`;
        if (r) html += `${dim('· MF')} <strong>${Math.round(r.mf)}</strong> ${dim('· RMF')} <strong>${Math.round(r.rmf)}</strong> ${dim(`· ${r.n_common}/${r.n_ref} íons da referência`)}`;
        msTitle.innerHTML = html;
        msTitle.title = msTitle.textContent;
        return;
    }

    if (msViewMode() === 'avg') {
        const a = s.avg;
        html = `<strong>Espectro médio</strong>${dim('·')} ${escapeHtml(s.name)}` +
            `${dim('·')} RT ${s.time[a.i0].toFixed(3)}–${s.time[a.i1].toFixed(3)} min` +
            dim(`· scans ${a.i0 + 1}–${a.i1 + 1} (${a.i1 - a.i0 + 1})`);
        if (a.bg && a.subtract) {
            html += dim(`· − fundo ${s.time[a.bg.i0].toFixed(3)}–${s.time[a.bg.i1].toFixed(3)} min`);
        }
    } else {
        const i = clampScan(s, s.scanIdx);
        html = `<strong>Espectro de massas</strong>${dim('·')} ${escapeHtml(s.name)}` +
            `${dim('·')} RT ${s.time[i].toFixed(4)} min` +
            dim(`· scan ${i + 1}/${s.time.length} · TIC ${formatBig(s.tic[i])}`);
    }

    if (spec && spec.mz.length) {
        let b = 0;
        for (let k = 1; k < spec.intensity.length; k++) if (spec.intensity[k] > spec.intensity[b]) b = k;
        html += `${dim('· pico base m/z')} ${+spec.mz[b].toFixed(2)}`;
    }
    msTitle.innerHTML = html;
    msTitle.title = msTitle.textContent;   // texto completo ao passar o mouse (o título é cortado quando falta espaço)
}

function showSpectrumMessage(text) {
    msEmptyMsg.style.display = text ? 'flex' : 'none';
    msEmptyMsg.innerText = text || '';
    if (text && msSpectrumDiv.data) { Plotly.purge(msSpectrumDiv); msSpectrumDiv._eicBound = false; }
}

async function loadSpectrum(s, idx) {
    idx = clampScan(s, idx);
    const token = ++msSpectrumToken;
    updateMsTitle();

    const key = `${s.id}|${idx}`;
    let spec = msSpectrumCache.get(key);
    if (!spec) {
        try {
            const res = await fetch(`/api/gcms/spectrum/${encodeURIComponent(s.id)}/${idx}`);
            if (!res.ok) {
                const err = await res.json().catch(() => ({}));
                throw new Error(err.detail || 'Falha ao obter o espectro.');
            }
            spec = await res.json();
            msSpectrumCache.set(key, spec);
            if (msSpectrumCache.size > 300) msSpectrumCache.delete(msSpectrumCache.keys().next().value);
        } catch (err) {
            if (token !== msSpectrumToken) return;
            msCurrentSpectrum = null;
            showSpectrumMessage(err.message);
            return;
        }
    }
    if (token !== msSpectrumToken) return;   // um scan mais novo já foi pedido

    msCurrentSpectrum = spec;
    showSpectrumMessage('');
    updateMsTitle();
    renderSpectrum();
}

// ── Rótulos de m/z do espectro ──
// Prioridade: o íon MAIS PESADO relevante (candidato a íon molecular) primeiro, depois os mais intensos.
// O espaçamento é medido em pixels (18 a 30 px por rótulo, conforme o modo), então só entra o que cabe sem sobrepor.
const LABEL_MODES = {
    few: { thr: 0.05, cap: 10, mplus: false, px: 30, text: 'poucos (≥ 5% do pico base)' },
    normal: { thr: 0.01, cap: 24, mplus: true, px: 26, text: 'normal (≥ 1%, mais o íon mais pesado)' },
    many: { thr: 0.002, cap: 60, mplus: true, px: 18, text: 'todos que couberem (≥ 0,2%, rótulos mais juntos)' }
};
let msLabelMode = 'normal';

function pickSpectrumLabels(mz, inten, maxI, xRange) {
    const cfg = LABEL_MODES[msLabelMode];
    const plotPx = Math.max(240, (msSpectrumDiv.clientWidth || 800) - 74);
    const sep = (xRange[1] - xRange[0]) * (cfg.px / plotPx);

    const cand = [];
    for (let i = 0; i < mz.length; i++) if (inten[i] >= maxI * cfg.thr) cand.push(i);
    cand.sort((a, b) => inten[b] - inten[a]);

    const order = [];
    if (cfg.mplus) {
        let heavy = -1;
        for (let i = 0; i < mz.length; i++) {
            if (inten[i] >= maxI * 0.003 && (heavy < 0 || mz[i] > mz[heavy])) heavy = i;
        }
        if (heavy >= 0) order.push(heavy);
    }
    cand.forEach(i => { if (!order.includes(i)) order.push(i); });

    const picked = [];
    for (const i of order) {
        if (picked.length >= cfg.cap) break;
        if (picked.every(j => Math.abs(mz[i] - mz[j]) >= sep)) picked.push(i);
    }
    return picked;
}

document.getElementById('ms-labels-btn').addEventListener('click', () => {
    const keys = Object.keys(LABEL_MODES);
    msLabelMode = keys[(keys.indexOf(msLabelMode) + 1) % keys.length];
    showToast(`Rótulos de m/z: ${LABEL_MODES[msLabelMode].text}`);
    renderSpectrum();
});

// Espectro em "sticks" (barras finas), com rótulo m/z nos picos mais relevantes
function renderSpectrum() {
    if (msPanelMode === 'cmp') { renderCompare(); return; }
    const spec = getDisplaySpectrum();
    const s = getActiveMs();
    if (!spec || !s || appMode !== 'gcms') return;

    if (spec.mz.length === 0) {
        showSpectrumMessage('Espectro vazio: a subtração do fundo zerou todos os íons.');
        return;
    }
    showSpectrumMessage('');

    const isDark = currentTheme === 'dark';
    const mz = spec.mz, inten = spec.intensity;
    let maxI = 0;
    for (const v of inten) if (v > maxI) maxI = v;
    if (!maxI) maxI = 1;
    const yv = msYRelative ? inten.map(v => (v / maxI) * 100) : inten;

    // Eixo m/z: do menor ao maior íon relevante (>= 0,5% do pico base) mais uma folga.
    // (Antes ficava fixo na faixa de aquisição do arquivo, deixando espaço vazio à direita.)
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < mz.length; i++) {
        if (inten[i] >= maxI * 0.005) {
            if (mz[i] < lo) lo = mz[i];
            if (mz[i] > hi) hi = mz[i];
        }
    }
    if (!isFinite(lo)) { lo = mz[0]; hi = mz[mz.length - 1]; }
    const span = hi - lo;
    const pad = Math.max(2, span * 0.04);
    const xRange = [Math.max(0, lo - pad), hi + pad];

    const sx = [], sy = [];
    for (let i = 0; i < mz.length; i++) { sx.push(mz[i], mz[i], null); sy.push(0, yv[i], null); }

    // Rótulos: ver pickSpectrumLabels (íon mais pesado primeiro, depois os mais intensos)
    const picked = pickSpectrumLabels(mz, inten, maxI, xRange);

    const traces = [
        { x: sx, y: sy, type: 'scatter', mode: 'lines', line: { color: s.color, width: 1.6 }, hoverinfo: 'none' },
        {
            x: mz, y: yv, type: 'scatter', mode: 'markers',
            marker: { size: 4, color: s.color },
            hovertemplate: `m/z %{x:.2f}<br>${msYRelative ? '%{y:.1f} %' : '%{y:,.2f}'}<br><i>clique: alternar EIC</i><extra></extra>`
        },
        {
            x: picked.map(i => mz[i]), y: picked.map(i => yv[i]),
            type: 'scatter', mode: 'text',
            text: picked.map(i => String(+mz[i].toFixed(1))),
            textposition: 'top center', cliponaxis: false,
            textfont: { size: 10, color: isDark ? '#f8fafc' : '#0f172a' },
            hoverinfo: 'none'
        }
    ];

    // Íons que estão no EIC: barras mais grossas na cor da curva correspondente
    msEic.ions.forEach((ion, k) => {
        const hx = [], hy = [];
        for (let i = 0; i < mz.length; i++) {
            if (Math.abs(mz[i] - ion) <= msEic.tol && yv[i] > 0) { hx.push(mz[i], mz[i], null); hy.push(0, yv[i], null); }
        }
        if (hx.length) {
            traces.push({
                x: hx, y: hy, type: 'scatter', mode: 'lines',
                line: { color: EIC_COLORS[k % EIC_COLORS.length], width: 3.5 }, hoverinfo: 'none'
            });
        }
    });

    const grid = isDark ? '#1a243b' : '#f1f5f9';
    const layout = {
        paper_bgcolor: isDark ? '#131b2e' : '#ffffff',
        plot_bgcolor: isDark ? '#131b2e' : '#ffffff',
        font: { color: isDark ? '#cbd5e1' : '#334155', size: 11 },
        margin: { l: 58, r: 16, t: 8, b: 36 },
        showlegend: false,
        hovermode: 'closest',
        xaxis: {
            title: { text: 'm/z', standoff: 4 }, gridcolor: grid, zeroline: false,
            range: xRange,
            // o zoom feito pelo usuário persiste ao trocar de scan (uirevision); sem zoom, o eixo segue cada espectro
            uirevision: `${s.id}|${spec.kind}`
        },
        yaxis: {
            title: { text: msYRelative ? 'Abundância relativa (%)' : 'Abundância', standoff: 4 },
            gridcolor: grid, zerolinecolor: isDark ? '#233152' : '#e2e8f0',
            rangemode: 'tozero',
            range: msYRelative ? [0, 118] : undefined,
            uirevision: msYRelative ? `rel|${spec.kind}` : `abs|${spec.key}`
        }
    };

    Plotly.react(msSpectrumDiv, traces, layout, { responsive: true, displaylogo: false });

    // clique numa barra: adiciona/remove o íon do EIC (o listener é registrado uma única vez por gráfico)
    if (!msSpectrumDiv._eicBound) {
        msSpectrumDiv.on('plotly_click', (ev) => {
            const pt = (ev.points || []).find(p => p.curveNumber === 1);   // 1 = marcadores das barras
            if (pt) toggleEicIon(pt.x);
        });
        msSpectrumDiv._eicBound = true;
    }
}

msRelBtn.addEventListener('click', () => {
    msYRelative = !msYRelative;
    msRelBtn.classList.toggle('active', msYRelative);
    msRelBtn.innerText = msYRelative ? '%' : 'abs';
    renderSpectrum();
});

// ── Espectro médio (intervalo escolhido no TIC) ──────────────────────────────
// Espectro realmente exibido na gaveta: um scan, ou a média do intervalo
// (já com o fundo subtraído, se ativado; valores negativos viram zero).
function getDisplaySpectrum() {
    const s = getActiveMs();
    if (!s) return null;

    if (msViewMode() === 'avg') {
        const r = msAvgResult;
        if (!r) return null;
        const sub = !!(s.avg.subtract && r.background);
        const mz = [], net = [], raw = [], bg = [];
        for (let i = 0; i < r.mz.length; i++) {
            const b = r.background ? r.background[i] : 0;
            const v = sub ? Math.max(r.intensity[i] - b, 0) : r.intensity[i];
            if (v > 0) { mz.push(r.mz[i]); net.push(v); raw.push(r.intensity[i]); bg.push(b); }
        }
        const bgKey = s.avg.bg ? `${s.avg.bg.i0}-${s.avg.bg.i1}` : '-';
        return {
            kind: 'avg', mz, intensity: net, raw, bg, subtracted: sub,
            key: `avg|${s.avg.i0}|${s.avg.i1}|${bgKey}|${s.avg.bin}|${sub}`
        };
    }

    const sp = msCurrentSpectrum;
    if (!sp) return null;
    return { kind: 'scan', mz: sp.mz, intensity: sp.intensity, key: `scan|${sp.scan}`, rt_min: sp.rt_min, scan: sp.scan };
}

// Mostra/esconde abas, navegação e controles da média conforme o modo da gaveta
function updateMsPanelUI() {
    const s = getActiveMs();
    const hasAvg = !!(s && s.avg);
    const hasRef = !!msReference;
    if (msPanelMode === 'avg' && !hasAvg) msPanelMode = 'scan';
    if (msPanelMode === 'cmp' && !hasRef) msPanelMode = 'scan';

    msModeTabs.style.display = (hasAvg || hasRef) ? 'flex' : 'none';
    msModeTabs.querySelectorAll('.ms-mode-tab').forEach(b => {
        const m = b.dataset.msmode;
        b.style.display = (m === 'avg' && !hasAvg) || (m === 'cmp' && !hasRef) ? 'none' : '';
        b.classList.toggle('active', m === msPanelMode);
    });
    msNav.style.display = msViewMode() === 'scan' ? 'flex' : 'none';
    msAvgControls.style.display = msPanelMode === 'avg' ? 'flex' : 'none';

    if (msPanelMode === 'avg') {
        msBinSelect.value = String(s.avg.bin);
        msBgToggle.style.display = s.avg.bg ? 'inline-block' : 'none';
        msBgToggle.classList.toggle('active', s.avg.subtract);
    }
}

async function refreshMsAvg() {
    const s = getActiveMs();
    if (!s || !s.avg) return;
    const a = s.avg;
    const token = ++msSpectrumToken;   // invalida pedidos de scan/média anteriores
    updateMsTitle();

    try {
        const res = await fetch('/api/gcms/average', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                sample_id: s.id, scan_start: a.i0, scan_end: a.i1, bin_width: a.bin,
                bg_start: a.bg ? a.bg.i0 : null, bg_end: a.bg ? a.bg.i1 : null
            })
        });
        if (!res.ok) {
            const err = await res.json().catch(() => ({}));
            throw new Error(err.detail || 'Falha ao calcular o espectro médio.');
        }
        const data = await res.json();
        if (token !== msSpectrumToken) return;
        msAvgResult = data;
    } catch (err) {
        if (token !== msSpectrumToken) return;
        msAvgResult = null;
        showSpectrumMessage(err.message);
        return;
    }
    showSpectrumMessage('');
    updateMsTitle();
    renderSpectrum();
}

// Chamado ao terminar a seleção início/fim no cromatograma
function applyMsInterval(purpose, t0, t1) {
    const s = getActiveMs();
    if (!s) return;
    const i0 = findClosestTimeIndex(s.time, t0);
    const i1 = findClosestTimeIndex(s.time, t1);

    if (purpose === 'ms-avg') {
        s.avg = {
            i0, i1,
            bg: s.avg ? s.avg.bg : null,
            bin: s.avg ? s.avg.bin : 1,
            subtract: s.avg ? s.avg.subtract : true
        };
    } else if (s.avg) {
        s.avg.bg = { i0, i1 };
        s.avg.subtract = true;
    }
    msPanelMode = 'avg';
    msAvgResult = null;
    updateApp();             // desenha as faixas (intervalo e fundo) no cromatograma
    showMsPanelIfNeeded();   // abre/atualiza a gaveta no modo Média
}

function toggleMsAvgMode() {
    if (appMode !== 'gcms') return;
    if (!getActiveMs()) { showToast('Importe um arquivo GC-MS primeiro.'); return; }
    togglePicking('ms-avg');
}

function toggleMsBgMode() {
    const s = getActiveMs();
    if (appMode !== 'gcms' || !s || !s.avg) return;
    togglePicking('ms-bg');
}

msAvgBtn.addEventListener('click', toggleMsAvgMode);
msBgBtn.addEventListener('click', toggleMsBgMode);

msModeTabs.querySelectorAll('.ms-mode-tab').forEach(btn => btn.addEventListener('click', () => {
    const s = getActiveMs();
    const mode = btn.dataset.msmode;
    if (!s || mode === msPanelMode || (mode === 'avg' && !s.avg) || (mode === 'cmp' && !msReference)) return;
    if (mode === 'cmp') { openCompareTab(); return; }
    msPanelMode = mode;
    updateMsPanelUI();
    updateApp();
    if (mode === 'avg') refreshMsAvg(); else loadSpectrum(s, s.scanIdx);
}));

msBinSelect.addEventListener('change', () => {
    const s = getActiveMs();
    if (!s || !s.avg) return;
    s.avg.bin = parseFloat(msBinSelect.value);
    msAvgResult = null;
    refreshMsAvg();
});

msBgToggle.addEventListener('click', () => {
    const s = getActiveMs();
    if (!s || !s.avg || !s.avg.bg) return;
    s.avg.subtract = !s.avg.subtract;
    updateMsPanelUI();
    updateMsTitle();
    renderSpectrum();
});

// ── Cromatograma de íon extraído (EIC) ──────────────────────────────────────
// Janela pela barra esquerda ou clique numa barra do espectro. O estado (íons, tolerância, TIC ao fundo)
// vale para todas as amostras: ao trocar de corrida o mesmo íon é recalculado, facilitando a comparação.
const eicBtn = document.getElementById('ms-eic-btn');
const eicModal = document.getElementById('eic-modal');
const eicIonsInput = document.getElementById('eic-ions');
const eicTolInput = document.getElementById('eic-tol');
const eicShowTicChk = document.getElementById('eic-show-tic');
let msEicToken = 0;

makeDraggable(eicModal, document.getElementById('eic-drag-handle'));

function syncEicInputs() {
    eicIonsInput.value = msEic.ions.map(x => +x.toFixed(2)).join(', ');
    eicTolInput.value = msEic.tol;
    eicShowTicChk.checked = msEic.showTic;
    eicBtn.classList.toggle('active-mode', msEic.ions.length > 0);
}

function closeEicModal() { eicModal.style.display = 'none'; }

eicBtn.addEventListener('click', () => {
    if (appMode !== 'gcms') return;
    if (eicModal.style.display !== 'none') { closeEicModal(); return; }
    syncEicInputs();
    eicModal.style.display = 'flex';
    eicIonsInput.focus();
});

// "43, 57; 91" -> [43, 57, 91] (sem repetidos, só valores válidos)
function parseIonList(text) {
    const out = [];
    text.split(/[\s,;]+/).forEach(tok => {
        const v = parseFloat(tok.replace(',', '.'));
        if (Number.isFinite(v) && v > 0 && v < 5000 && !out.some(x => Math.abs(x - v) < 1e-6)) out.push(+v.toFixed(4));
    });
    return out;
}

async function refreshMsEic() {
    const s = getActiveMs();
    if (!s || msEic.ions.length === 0) {
        msEicData = null;
        updateApp();
        renderBottomPanel();
        renderSpectrum();
        return;
    }
    if (isEicActive()) { updateApp(); renderBottomPanel(); renderSpectrum(); return; }

    const token = ++msEicToken;
    const key = eicKey(s);
    let data;
    try {
        const res = await fetch('/api/gcms/eic', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ sample_id: s.id, ions: msEic.ions, tol: msEic.tol })
        });
        if (!res.ok) {
            const err = await res.json().catch(() => ({}));
            throw new Error(err.detail || 'Falha ao calcular o EIC.');
        }
        data = await res.json();
    } catch (err) {
        if (token === msEicToken) { msEicData = null; showToast(err.message); }
        return;
    }
    if (token !== msEicToken) return;   // um pedido mais novo (ou outra amostra) já substituiu este
    msEicData = {
        key, series: data.series,
        signals: data.series.map((arr, k) => makeMsSignal(s, msEic.ions[k], msEic.tol, EIC_COLORS[k % EIC_COLORS.length], arr))
    };
    updateApp();
    renderBottomPanel();
    msAutorange();
    renderSpectrum();
}

async function applyEicFromModal() {
    const ions = parseIonList(eicIonsInput.value);
    if (ions.length === 0) { showToast('Informe ao menos um m/z (ex.: 43, 57).'); return; }
    if (ions.length > EIC_MAX) { ions.length = EIC_MAX; showToast(`Máximo de ${EIC_MAX} íons: os demais foram ignorados.`); }
    let tol = parseFloat(eicTolInput.value);
    if (!(tol >= 0.01)) tol = 0.5;
    msEic = { ions, tol: Math.min(5, tol), showTic: eicShowTicChk.checked };
    syncEicInputs();
    await refreshMsEic();
}

function clearEic() {
    msEic.ions = [];
    msEicData = null;
    msEicToken++;
    syncEicInputs();
    updateApp();
    renderBottomPanel();
    msAutorange();
    renderSpectrum();
}

// Clique numa barra do espectro: adiciona o íon ao EIC ou, se já estiver, remove
function toggleEicIon(mz) {
    if (appMode !== 'gcms' || !getActiveMs()) return;
    const ion = msEic.tol >= 0.5 ? Math.round(mz) : +mz.toFixed(2);
    const idx = msEic.ions.findIndex(x => Math.abs(x - ion) < 1e-6);
    if (idx >= 0) {
        msEic.ions.splice(idx, 1);
        showToast(`EIC: m/z ${ion} removido`);
    } else {
        if (msEic.ions.length >= EIC_MAX) { showToast(`Máximo de ${EIC_MAX} íons no EIC.`); return; }
        msEic.ions.push(ion);
        showToast(`EIC: m/z ${ion} adicionado`);
    }
    msEicData = null;
    syncEicInputs();
    refreshMsEic().then(() => { if (msEic.ions.length === 0) msAutorange(); });
}

document.getElementById('btn-eic-apply').addEventListener('click', applyEicFromModal);
document.getElementById('btn-eic-clear').addEventListener('click', clearEic);
document.getElementById('eic-close-btn').addEventListener('click', closeEicModal);
[eicIonsInput, eicTolInput].forEach(inp => inp.addEventListener('keydown', (e) => { if (e.key === 'Enter') applyEicFromModal(); }));
eicShowTicChk.addEventListener('change', () => { msEic.showTic = eicShowTicChk.checked; updateApp(); });

// ── Comparação de espectros (aba "Comparar" + janela ⇆ Ref.) ─────────────────
// Compara o espectro exibido (um scan ou a média) com uma referência: outro espectro seu ou uma entrada
// de uma biblioteca MSP. As notas (MF direta e RMF reversa, 0–999) vêm do servidor.
const specrefBtn = document.getElementById('ms-ref-btn');
const specrefModal = document.getElementById('specref-modal');
const refInfoEl = document.getElementById('ref-info');
const libSelect = document.getElementById('lib-select');
const libFile = document.getElementById('lib-file');
const libResults = document.getElementById('lib-results');
const libNameQ = document.getElementById('lib-name-q');
let msCmpToken = 0;

makeDraggable(specrefModal, document.getElementById('specref-drag-handle'));

// Qual espectro está em exibição: o da aba atual (na aba Comparar, o escolhido em msCmpSource)
function msViewMode() {
    const s = getActiveMs();
    const m = msPanelMode === 'cmp' ? msCmpSource : msPanelMode;
    return (m === 'avg' && s && s.avg) ? 'avg' : 'scan';
}

async function renderCompare() {
    const spec = getDisplaySpectrum();
    const s = getActiveMs();
    if (!spec || !s || !msReference || appMode !== 'gcms') return;
    if (spec.mz.length === 0) { showSpectrumMessage('Espectro vazio: nada a comparar.'); return; }

    const token = ++msCmpToken;
    let res;
    try {
        const r = await fetch('/api/gcms/compare', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ mz: spec.mz, intensity: spec.intensity, ref_mz: msReference.mz, ref_intensity: msReference.intensity })
        });
        if (!r.ok) throw new Error('Falha ao comparar os espectros.');
        res = await r.json();
    } catch (err) {
        if (token === msCmpToken) showSpectrumMessage(err.message);
        return;
    }
    if (token !== msCmpToken || msPanelMode !== 'cmp') return;
    showSpectrumMessage('');
    msCmpResult = res;
    updateMsTitle();

    const isDark = currentTheme === 'dark';
    const refColor = '#f59e0b';
    const mz = res.mz, up = res.query, dn = res.ref;
    const ux = [], uy = [], dx = [], dy = [];
    for (let i = 0; i < mz.length; i++) {
        if (up[i] > 0) { ux.push(mz[i], mz[i], null); uy.push(0, up[i], null); }
        if (dn[i] > 0) { dx.push(mz[i], mz[i], null); dy.push(0, -dn[i], null); }
    }

    // rótulos (mesmo critério do espectro simples; ambos os lados já estão normalizados a 100)
    const sigX = [];
    for (let i = 0; i < mz.length; i++) if (up[i] >= 0.5 || dn[i] >= 0.5) sigX.push(mz[i]);
    const lo0 = sigX.length ? sigX[0] : mz[0], hi0 = sigX.length ? sigX[sigX.length - 1] : mz[mz.length - 1];
    const pad0 = Math.max(2, (hi0 - lo0) * 0.04);
    const xr = [Math.max(0, lo0 - pad0), hi0 + pad0];
    const labels = (arr, sign) => {
        const picked = pickSpectrumLabels(mz, arr, 100, xr);
        return { x: picked.map(i => mz[i]), y: picked.map(i => sign * arr[i]), text: picked.map(i => String(+mz[i].toFixed(0))) };
    };
    const lu = labels(up, 1), ld = labels(dn, -1);
    const txtColor = isDark ? '#f8fafc' : '#0f172a';
    const traces = [
        { x: ux, y: uy, type: 'scatter', mode: 'lines', line: { color: s.color, width: 2 }, hoverinfo: 'none' },
        { x: dx, y: dy, type: 'scatter', mode: 'lines', line: { color: refColor, width: 2 }, hoverinfo: 'none' },
        { x: mz.filter((_, i) => up[i] > 0), y: up.filter(v => v > 0), type: 'scatter', mode: 'markers', marker: { size: 4, color: s.color }, hovertemplate: 'm/z %{x}<br>%{y:.1f} %<extra>amostra</extra>' },
        { x: mz.filter((_, i) => dn[i] > 0), y: dn.filter(v => v > 0).map(v => -v), type: 'scatter', mode: 'markers', marker: { size: 4, color: refColor }, hovertemplate: 'm/z %{x}<br>%{customdata:.1f} %<extra>referência</extra>', customdata: dn.filter(v => v > 0) },
        { x: lu.x, y: lu.y, type: 'scatter', mode: 'text', text: lu.text, textposition: 'top center', textfont: { size: 10, color: txtColor }, hoverinfo: 'none', cliponaxis: false },
        { x: ld.x, y: ld.y, type: 'scatter', mode: 'text', text: ld.text, textposition: 'bottom center', textfont: { size: 10, color: txtColor }, hoverinfo: 'none', cliponaxis: false }
    ];

    const sig = [];
    for (let i = 0; i < mz.length; i++) if (up[i] >= 0.5 || dn[i] >= 0.5) sig.push(mz[i]);
    const lo = sig.length ? sig[0] : mz[0], hi = sig.length ? sig[sig.length - 1] : mz[mz.length - 1];
    const pad = Math.max(2, (hi - lo) * 0.04);
    const grid = isDark ? '#1a243b' : '#f1f5f9';
    const layout = {
        paper_bgcolor: isDark ? '#131b2e' : '#ffffff',
        plot_bgcolor: isDark ? '#131b2e' : '#ffffff',
        font: { color: isDark ? '#cbd5e1' : '#334155', size: 11 },
        margin: { l: 58, r: 16, t: 8, b: 36 },
        showlegend: false, hovermode: 'closest',
        xaxis: { title: { text: 'm/z', standoff: 4 }, gridcolor: grid, zeroline: false, range: [Math.max(0, lo - pad), hi + pad], uirevision: `${s.id}|cmp` },
        yaxis: {
            title: { text: 'Abundância relativa (%)', standoff: 4 }, gridcolor: grid, range: [-122, 122],
            zeroline: true, zerolinecolor: isDark ? '#64748b' : '#94a3b8', zerolinewidth: 1.5,
            tickvals: [-100, -50, 0, 50, 100], ticktext: ['100', '50', '0', '50', '100'], uirevision: 'cmp-y'
        },
        annotations: [
            {
                xref: 'paper', yref: 'y', x: 0.005, y: 112, xanchor: 'left', showarrow: false, font: { size: 11, color: txtColor },
                text: `<b>MF ${Math.round(res.mf)}</b> · RMF ${Math.round(res.rmf)} · ${res.n_common}/${res.n_ref} íons da referência`
            },
            { xref: 'paper', yref: 'y', x: 0.995, y: 112, text: escapeHtml(s.name), showarrow: false, xanchor: 'right', font: { size: 10, color: s.color } },
            { xref: 'paper', yref: 'y', x: 0.995, y: -112, text: escapeHtml(msReference.name), showarrow: false, xanchor: 'right', font: { size: 10, color: refColor } }
        ]
    };
    Plotly.react(msSpectrumDiv, traces, layout, { responsive: true, displaylogo: false });
    msSpectrumDiv._eicBound = msSpectrumDiv._eicBound || false;
}

// ── Referência e biblioteca ──
function updateRefInfo() {
    refInfoEl.textContent = msReference
        ? `${msReference.name}  (${msReference.mz.length} íons · ${msReference.source})`
        : 'Nenhuma referência';
    refInfoEl.title = refInfoEl.textContent;
}

function setReference(ref) {
    msReference = ref;
    msCmpResult = null;
    updateRefInfo();
    updateMsPanelUI();
}

function openCompareTab() {
    const s = getActiveMs();
    if (!s || !msReference) return;
    if (msPanelMode !== 'cmp') msCmpSource = (msPanelMode === 'avg' && s.avg) ? 'avg' : 'scan';
    msPanelMode = 'cmp';
    updateMsPanelUI();
    updateApp();
    if (msCmpSource === 'avg') refreshMsAvg(); else loadSpectrum(s, s.scanIdx);
}

function leaveCompareIfNoRef() {
    if (!msReference && msPanelMode === 'cmp') {
        msPanelMode = 'scan';
        updateMsPanelUI();
        updateApp();
        const s = getActiveMs();
        if (s) loadSpectrum(s, s.scanIdx);
    }
}

async function refreshLibraries() {
    try {
        const res = await fetch('/api/gcms/library/list');
        const data = await res.json();
        const keep = libSelect.value;
        libSelect.innerHTML = data.libraries.length
            ? `<option value="">Todas (${data.libraries.length})</option>` + data.libraries.map(l => `<option value="${escapeHtml(l.lib_id)}">${escapeHtml(l.name)} — ${l.n_entries} espectros</option>`).join('')
            : '<option value="">(carregue um arquivo .msp)</option>';
        if ([...libSelect.options].some(o => o.value === keep)) libSelect.value = keep;
    } catch (err) { /* servidor indisponível: mantém a lista anterior */ }
}

function renderLibHits(hits, withScores) {
    libResults.innerHTML = '';
    if (!hits.length) {
        libResults.innerHTML = '<tr><td colspan="5" class="muted">Nada encontrado.</td></tr>';
        return;
    }
    hits.forEach((h, i) => {
        const tr = document.createElement('tr');
        tr.innerHTML = `<td>${i + 1}</td><td title="${escapeHtml(h.name)}">${escapeHtml(h.name.length > 34 ? h.name.slice(0, 33) + '…' : h.name)}</td>
            <td class="score">${withScores ? Math.round(h.mf) : '—'}</td><td class="score">${withScores ? Math.round(h.rmf) : '—'}</td>
            <td><button class="ms-bg-btn">Usar</button></td>`;
        tr.querySelector('button').addEventListener('click', () => useLibraryHit(h));
        libResults.appendChild(tr);
    });
}

async function useLibraryHit(h) {
    try {
        const res = await fetch(`/api/gcms/library/entry/${encodeURIComponent(h.lib_id)}/${h.index}`);
        if (!res.ok) throw new Error('Entrada não encontrada.');
        const e = await res.json();
        setReference({ name: e.name, mz: e.mz, intensity: e.intensity, source: `biblioteca ${h.lib_id}`, meta: e.meta || {} });
        openCompareTab();
        showToast(`Referência: ${e.name}`);
    } catch (err) {
        showToast(err.message);
    }
}

document.getElementById('lib-upload-btn').addEventListener('click', () => libFile.click());
libFile.addEventListener('change', async () => {
    if (!libFile.files.length) return;
    const fd = new FormData();
    for (const f of libFile.files) fd.append('files', f);
    try {
        const res = await fetch('/api/gcms/library/upload', { method: 'POST', body: fd });
        const data = await res.json();
        if (data.errors && data.errors.length) alert('Não foi possível ler:\n' + data.errors.map(e => `• ${e.filename}: ${e.error}`).join('\n'));
        if (data.success && data.success.length) {
            await refreshLibraries();
            libSelect.value = data.success[0].lib_id;
            showToast(data.success.map(l => `${l.name}: ${l.n_entries} espectros`).join(' · '));
        }
    } catch (err) { alert('Falha ao enviar a biblioteca.'); }
    libFile.value = '';
});

document.getElementById('lib-search').addEventListener('click', async () => {
    const spec = getDisplaySpectrum();
    if (!spec || !spec.mz.length) { showToast('Não há espectro exibido para buscar.'); return; }
    try {
        const res = await fetch('/api/gcms/library/search', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ lib_id: libSelect.value || null, mz: spec.mz, intensity: spec.intensity, top: 10 })
        });
        if (!res.ok) throw new Error('Carregue uma biblioteca .msp primeiro.');
        renderLibHits((await res.json()).hits, true);
    } catch (err) { showToast(err.message); }
});

async function findLibraryByName() {
    if (!libNameQ.value.trim()) return;
    try {
        const res = await fetch('/api/gcms/library/find', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ lib_id: libSelect.value || null, q: libNameQ.value })
        });
        renderLibHits((await res.json()).hits, false);
    } catch (err) { showToast('Falha na busca.'); }
}
document.getElementById('lib-find').addEventListener('click', findLibraryByName);
libNameQ.addEventListener('keydown', (e) => { if (e.key === 'Enter') findLibraryByName(); });

document.getElementById('ref-set-current').addEventListener('click', () => {
    const s = getActiveMs();
    const spec = getDisplaySpectrum();
    if (!s || !spec || !spec.mz.length) { showToast('Não há espectro exibido para fixar.'); return; }
    const label = spec.kind === 'avg'
        ? `${s.name} · média ${s.time[s.avg.i0].toFixed(3)}–${s.time[s.avg.i1].toFixed(3)} min`
        : `${s.name} · RT ${spec.rt_min.toFixed(3)} min`;
    setReference({ name: label, mz: spec.mz.slice(), intensity: spec.intensity.slice(), source: 'amostra', meta: {} });
    showToast('Referência fixada. Vá para outro scan/amostra e use a aba Comparar.');
});

document.getElementById('ref-clear').addEventListener('click', () => {
    setReference(null);
    leaveCompareIfNoRef();
    libResults.innerHTML = '';
});

function closeSpecrefModal() { specrefModal.style.display = 'none'; }
specrefBtn.addEventListener('click', () => {
    if (specrefModal.style.display !== 'none') { closeSpecrefModal(); return; }
    if (appMode !== 'gcms') return;
    updateRefInfo();
    refreshLibraries();
    specrefModal.style.display = 'flex';
});
document.getElementById('specref-close-btn').addEventListener('click', closeSpecrefModal);

// ── Exportação (CSV e MSP) do espectro exibido ──
function downloadText(text, filename, mime) {
    const blob = new Blob([text], { type: mime });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}

// RT de referência do espectro exibido: o scan, ou o topo do TIC dentro do intervalo médio
function getMsExportInfo(s, spec) {
    const safe = s.name.replace(/\s+/g, '_');
    if (spec.kind === 'avg') {
        const a = s.avg;
        let apex = a.i0;
        for (let i = a.i0; i <= a.i1; i++) if (s.tic[i] > s.tic[apex]) apex = i;
        return { rt: s.time[apex], safe, label: `${safe}_RT${s.time[apex].toFixed(3)}_media` };
    }
    return { rt: spec.rt_min, safe, label: `${safe}_RT${spec.rt_min.toFixed(3)}` };
}

msExportBtn.addEventListener('click', () => {
    const s = getActiveMs();
    const spec = getDisplaySpectrum();
    if (!s || !spec || !spec.mz.length) { alert('Nenhum espectro para exportar.'); return; }
    const info = getMsExportInfo(s, spec);
    let csv = '\uFEFF';
    if (spec.kind === 'avg') {
        const a = s.avg;
        csv += `# ${s.name};espectro médio;scans ${a.i0 + 1}-${a.i1 + 1};RT ${s.time[a.i0]}-${s.time[a.i1]} min;bin ${a.bin}` +
            (a.bg ? `;fundo scans ${a.bg.i0 + 1}-${a.bg.i1 + 1} ${spec.subtracted ? '(subtraído)' : '(não subtraído)'}` : '') + '\n';
        csv += 'm/z;Intensidade;Média bruta;Fundo\n';
        for (let i = 0; i < spec.mz.length; i++) csv += `${spec.mz[i]};${spec.intensity[i]};${spec.raw[i]};${spec.bg[i]}\n`;
    } else {
        csv += `# ${s.name};RT (min);${spec.rt_min};scan;${spec.scan + 1}\n` + 'm/z;Intensidade\n';
        for (let i = 0; i < spec.mz.length; i++) csv += `${spec.mz[i]};${spec.intensity[i]}\n`;
    }
    downloadText(csv, `espectro_${info.label}.csv`, 'text/csv;charset=utf-8;');
});

// MSP: formato de texto lido por bibliotecas/buscas espectrais (NIST MS Search, MS-DIAL, etc.).
// Intensidades normalizadas para 999 (convenção NIST); íons que arredondam para 0 são omitidos.
msExportMspBtn.addEventListener('click', () => {
    const s = getActiveMs();
    const spec = getDisplaySpectrum();
    if (!s || !spec || !spec.mz.length) { alert('Nenhum espectro para exportar.'); return; }

    let maxI = 0;
    for (const v of spec.intensity) if (v > maxI) maxI = v;
    const peaks = [];
    for (let i = 0; i < spec.mz.length; i++) {
        const v = Math.round((spec.intensity[i] / maxI) * 999);
        if (v > 0) peaks.push(`${Number.isInteger(spec.mz[i]) ? spec.mz[i] : +spec.mz[i].toFixed(4)} ${v};`);
    }
    if (peaks.length === 0) { alert('Espectro vazio.'); return; }

    const info = getMsExportInfo(s, spec);
    const ascii = (t) => String(t).normalize('NFD').replace(/[\u0300-\u036f]/g, '').replace(/[\r\n]+/g, ' ');
    let comment;
    if (spec.kind === 'avg') {
        const a = s.avg;
        comment = `Espectro medio dos scans ${a.i0 + 1}-${a.i1 + 1} (n=${a.i1 - a.i0 + 1}), RT ${s.time[a.i0].toFixed(3)}-${s.time[a.i1].toFixed(3)} min; ` +
            `m/z alinhado em bins de ${a.bin}; ` +
            (spec.subtracted ? `fundo subtraido (scans ${a.bg.i0 + 1}-${a.bg.i1 + 1}); ` : 'sem subtracao de fundo; ');
    } else {
        comment = `Scan ${spec.scan + 1}; `;
    }
    comment += `amostra ${ascii(s.name)}; RT em minutos; intensidades normalizadas para 999`;

    const name = ascii(spec.kind === 'avg' ? `${s.name} RT ${info.rt.toFixed(3)} (media)` : `${s.name} RT ${info.rt.toFixed(3)}`);
    const lines = [];
    for (let i = 0; i < peaks.length; i += 5) lines.push(peaks.slice(i, i + 5).join(' '));

    const msp = `Name: ${name}\nComments: ${comment}\nRT: ${info.rt.toFixed(4)}\nNum Peaks: ${peaks.length}\n${lines.join('\n')}\n\n`;
    downloadText(msp, `${info.label}.msp`, 'text/plain;charset=utf-8;');
});


msMinimizeBtn.addEventListener('click', () => {
    const isMin = msPanel.classList.toggle('minimized');
    msMinimizeBtn.innerText = isMin ? '□' : '_';
    setTimeout(() => { Plotly.Plots.resize(plotDiv); if (msSpectrumDiv.data) Plotly.Plots.resize(msSpectrumDiv); }, 30);
});

// Redimensionamento vertical da gaveta do espectro
let isResizingMs = false;
msPanelResizer.addEventListener('mousedown', () => {
    isResizingMs = true;
    msPanelResizer.classList.add('resizing');
    document.body.style.cursor = 'row-resize';
    document.body.style.userSelect = 'none';
});
document.addEventListener('mousemove', (e) => {
    if (!isResizingMs) return;
    const newHeight = window.innerHeight - e.clientY;
    if (newHeight >= 110 && newHeight <= window.innerHeight * 0.75) {
        msPanel.classList.remove('minimized');
        msMinimizeBtn.innerText = '_';
        msPanel.style.height = `${newHeight}px`;
        Plotly.Plots.resize(plotDiv);
    }
});
document.addEventListener('mouseup', () => {
    if (isResizingMs) {
        isResizingMs = false;
        msPanelResizer.classList.remove('resizing');
        document.body.style.cursor = '';
        document.body.style.userSelect = '';
        Plotly.Plots.resize(plotDiv);
    }
});

// O espectro acompanha qualquer mudança de tamanho do seu container (gaveta, barras laterais, janela)
new ResizeObserver(() => { if (msSpectrumDiv.data) Plotly.Plots.resize(msSpectrumDiv); }).observe(msSpectrumDiv);

// Inicialização
fetchProjectsInfo();
updateCurrentProjectDisplay();
renderRampModal();