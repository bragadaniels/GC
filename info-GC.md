# INFO.md — Guia de Configuração e Uso do GCAnalyzer

Este documento detalha todos os parâmetros ajustáveis da classe `ProcessingMethod` no arquivo `gc.py` e fornece exemplos de como utilizar as classes principais (`GCAnalyzer`, `GCReport`, `GCVisualizer`).

---

## 🟢 1. Parâmetros de Processamento (`ProcessingMethod`)

A classe `ProcessingMethod` controla todo o pipeline. Os parâmetros estão agrupados por função.

### A. Linha de Base (Baseline) e Ruído

Responsável por separar o sinal químico do "fundo" eletrônico ou do gradiente do solvente.

| Parâmetro | Padrão | Descrição e Comportamento |
| :--- | :--- | :--- |
| **`baseline_lam`** | `1e7` | **Suavidade da Baseline (Lambda).**<br>• **Aumentar ($>10^9$):** Baseline fica rígida (reta). *Risco:* Cortar picos largos ou ignorar subida de gradiente.<br>• **Diminuir ($<10^5$):** Baseline fica flexível. *Risco:* Entrar dentro dos picos e subtrair área do analito. |
| **`baseline_p`** | `0.0001` | **Assimetria.** Peso dado aos pontos acima da baseline.<br>• **Aumentar:** A baseline sobe. *Risco:* Cortar a base dos picos.<br>• **Diminuir:** A baseline desce. *Risco:* Ficar muito abaixo do ruído, gerando áreas excessivas.<br>⚡ **Acoplamento:** Se diminuir `lam` (mais flexível), deve diminuir `p` para evitar que a baseline invada os picos. |
| **`noise_percentile`** | `20` | **Estimativa de Ruído.** Percentil da intensidade usado para calcular o ruído de fundo.<br>• **Aumentar:** Em cromatogramas muito cheios, pode considerar picos pequenos como ruído.<br>• **Ideal:** 10 a 20 para a maioria dos casos. |

### B. Pré-Processamento e Filtros Temporais

Definem *onde* e *o que* o algoritmo deve olhar antes mesmo de tentar detectar picos.

| Parâmetro | Padrão | Descrição e Comportamento |
| :--- | :--- | :--- |
| **`t_start_integration`** | `0.0` | **Corte Inicial (segundos).** Ignora completamente qualquer sinal antes deste tempo.<br>• **Uso:** Remover o pico gigante do solvente ou instabilidade inicial do detector que quebra a escala.<br>• **Ex:** `60.0` (ignora o primeiro minuto). |
| **`solvent_rt_cutoff_s`** | `60.0` | **Filtro Pós-Detecção.** Picos detectados antes deste tempo são descartados.<br>• **Diferença para t_start:** `t_start` corta o sinal bruto; este parâmetro corta picos já integrados. |
| **`solvent_area_factor`** | `5.0` | **Filtro de Solvente por Área.** Se `Area_Pico > Fator * Mediana_Areas`, o pico é removido (exceto se for Padrão Interno).<br>• **Uso:** Remover "bolhas" ou picos gigantes de retrolavagem. |
| **`min_area_threshold`** | `0.0` | **Área Mínima Absoluta.**<br>• **Uso:** Descartar "spikes" (ruído elétrico) que passam no filtro de altura (SNR) mas não têm massa (área) real.<br>• **Ex:** `100.0`. |

### C. Detecção de Picos e Suavização

Parâmetros que determinam se uma oscilação no sinal é um pico ou não.

| Parâmetro | Padrão | Descrição e Comportamento |
| :--- | :--- | :--- |
| **`snr_threshold`** | `3.0` | **Sinal-Ruído Mínimo.**<br>• **Aumentar:** Mais rigoroso. Perde picos traço.<br>• **Diminuir:** Detecta mais picos. *Risco:* Integrar ruído como pico. |
| **`min_width_seconds`** | `1.0` | **Largura Mínima (s).** Picos mais estreitos que isso são ignorados.<br>• **Uso:** Filtrar ruído de alta frequência. |
| **`min_distance_seconds`** | `2.0` | **Distância Mínima entre Ápices.**<br>• **Uso:** Evita detectar múltiplos topos em um pico largo e ruidoso (serrilhado). |
| **`sg_window_length`** | `0` (Auto) | **Suavização Savitzky-Golay (Janela).** Número de pontos para suavizar (deve ser ímpar).<br>• **0:** Calcula automático (~1.5x a largura média).<br>• **Aumentar:** Cromatograma mais liso, SNR melhora. *Risco:* Achatar picos e perder resolução (Rs). |
| **`sg_polyorder`** | `4` | **Ordem do Polinômio SG.**<br>• **Ideal:** 3 ou 4.<br>⚡ **Acoplamento:** Deve ser sempre menor que `sg_window_length`. |
| **`slope_threshold_factor`** | `0.10` | **Sensibilidade da Inclinação (Início/Fim).** Define onde a integração começa/termina baseada na derivada.<br>• **Aumentar (0.5):** Integração termina cedo. Corta cauda (tailing).<br>• **Diminuir (0.01):** Integração vai longe na baseline. *Risco:* Fundir picos vizinhos. |

### D. Controle Manual e Avançado (Forçar/Inibir)

Parâmetros para quando a detecção automática falha ou precisa de intervenção humana.

| Parâmetro | Padrão | Descrição e Comportamento |
| :--- | :--- | :--- |
| **`integration_inhibit_windows`** | `[]` | **Zonas Mortas.** Lista de tuplas `[(inicio, fim)]`.<br>• **O que faz:** Impede a detecção de *qualquer* pico automático nestas janelas.<br>• **Uso:** Mascarar troca de válvulas ou artefatos conhecidos.<br>• **Ex:** `[(120.0, 125.5), (300, 310)]`. |
| **`force_integration_windows`** | `[]` | **Integração Forçada.** Lista de tuplas `[(inicio, fim)]`.<br>• **O que faz:** Força a criação de um pico nesta janela, ignorando SNR, limiar de área ou filtros de solvente.<br>• **Como funciona:** Acha o máximo local na janela e ajusta um modelo.<br>• **Uso:** Quantificar impurezas conhecidas que estão abaixo do limite de detecção automático.<br>⚡ **Atenção:** Bypass total de filtros de segurança. |
| **`expected_peaks_count`** | `None` | **QC de Contagem.**<br>• **O que faz:** Apenas gera um aviso (WARN) no log se o número de picos encontrados for diferente deste valor. Não altera o processamento. |

### E. Lógica de Separação de Picos (Coeluição)

Árvore de decisão para picos sobrepostos (quando $Rs <$ `rs_deconv_threshold`).

| Parâmetro | Padrão | Descrição |
| :--- | :--- | :--- |
| **`rs_deconv_threshold`** | `1.2` | Se Resolução > isso, picos são tratados como isolados. Se menor, analisa o vale. |
| **`valley_pct_independent`** | `85.0` | Se o vale desce até 85% da altura, separa com corte vertical (Baseline-to-Baseline). |
| **`valley_pct_dropline`** | `50.0` | Se o vale é profundo (>50%), usa **Drop-Line** (linha reta do vale à baseline). |
| **`valley_pct_skim_max`** | `25.0` | Se o vale é raso (<25%) e a diferença de altura é grande, tenta **Tangent Skim**. |
| **`height_ratio_rider`** | `0.15` | Para Skim: o pico menor ("rider") deve ser no máximo 15% do maior ("parent"). Caso contrário, usa Deconvolução. |

---

## 💻 2. Exemplos de Uso das Classes Principais

### Exemplo 1: `GCAnalyzer` (O Motor de Processamento)

Use para processar arquivos `.cdf` brutos e obter DataFrames de resultados.

```python
from gc_pipeline import GCAnalyzer, ProcessingMethod

# 1. Configurar o método (pode carregar de JSON também)
metodo = ProcessingMethod(
    name="Metodo_Oleo_Essencial",
    snr_threshold=5.0,
    min_width_seconds=0.8,
    # Forçar a integração de um pico conhecido em 145s que é muito fraco
    force_integration_windows=[(144.5, 146.0)]
)

# 2. Instanciar o analisador
analyzer = GCAnalyzer(method=metodo, run_id="Lote_2023_A")

# 3. Processar um arquivo (Single Run)
try:
    # Retorna vetores (tempo, sinal) e o DataFrame de resultados
    rt, raw_signal = analyzer.read_cdf("dados/amostra_01.cdf")
    corrected_signal, baseline = analyzer.remove_baseline(rt, raw_signal)
    df_resultados = analyzer.integrate(rt, corrected_signal)
    
    # Calcular métricas farmacopeicas (USP/EP)
    df_completo = analyzer.compute_usp_metrics(rt, corrected_signal, df_resultados)
    
    print(df_completo[["rt", "area", "N_plates_ep", "tailing_factor_usp"]])

except Exception as e:
    print(f"Erro no processamento: {e}")

# 4. Exportar Audit Trail (Rastreabilidade)
analyzer.export_audit("audit_trail.json")
```

---

### Exemplo 2: `GCReport` (Geração de PDF)

Use para criar relatórios laudos profissionais em PDF.

```python
from gc_pipeline import GCReport

# Supondo que você já rodou o analyzer acima e tem os dados:
# rt, raw, corrected, baseline, df_completo

# 1. Criar o relatório
report = GCReport(
    analyzer=analyzer, # Passa o analyzer para pegar metadados do método
    title="Laudo de Análise de Pureza",
    analyst="Dr. Silva",
    sample_info="Lote XYZ-99, Vial 2"
)

# 2. Adicionar a corrida (pode adicionar várias para comparação)
report.add_run(
    rt=rt,
    raw=raw_signal, # Opcional, mas bom para histórico
    corrected=corrected_signal,
    baseline=baseline,
    results_df=df_completo,
    label="Injeção 01"
)

# 3. Gerar o PDF
pdf_path = report.build("Relatorio_Final.pdf")
print(f"PDF gerado em: {pdf_path}")
```

---

### Exemplo 3: `GCVisualizer` (Gráficos Interativos)

Use para análise exploratória em HTML/Plotly. Ótimo para "debugar" métodos.

```python
from gc_pipeline import GCVisualizer

# 1. Instanciar visualizador
viz = GCVisualizer(show_on_plot=True) # show_on_plot=True abre o navegador automaticamente

# 2. Adicionar dados (pode ser de várias corridas para sobreposição)
viz.add_run(
    label="Amostra A",
    rt=rt,
    signal=corrected_signal, # Sinal corrigido (sem baseline)
    results_df=df_completo
)

# 3. Gerar plots
# Plot Individual (com zoom, hover nos picos)
viz.plot_single("Amostra A", filename="grafico_amostra_a.html")

# Se tivesse adicionado mais corridas:
# viz.plot_overlay(title="Comparativo Lote A vs B", normalize=True)
# viz.plot_stacked(offset_factor=1.2)
```

---

# 🛠️ GCAnalyzer: Guia de Solução de Problemas (Cheat Sheet)

Use esta tabela para diagnosticar problemas visuais no seu cromatograma e ajustar o parâmetro correto no `ProcessingMethod`.

---

## 📉 1. Problemas de Linha de Base (Baseline)

O "chão" do seu cromatograma está errado.

| Sintoma Visual | Diagnóstico | Ação Recomendada | Parâmetro |
| :--- | :--- | :--- | :--- |
| **Baseline corta o meio de picos largos** | A baseline está muito "rígida" e não acompanha a subida do pico. | **Diminuir** a rigidez. | `baseline_lam` (ex: $10^7 \to 10^5$) |
| **Baseline "entra" na base dos picos** | A baseline está muito alta. | **Diminuir** a assimetria. | `baseline_p` (ex: $0.001 \to 0.0001$) |
| **Baseline flutua muito abaixo do sinal** | A baseline está sendo empurrada para baixo por ruídos negativos. | **Aumentar** a assimetria. | `baseline_p` (ex: $10^{-5} \to 10^{-3}$) |
| **Baseline segue o ruído (ondulada)** | A baseline está muito "mole" ou flexível. | **Aumentar** a rigidez. | `baseline_lam` (ex: $10^5 \to 10^8$) |

> **💡 Dica de Ouro:** Para cromatogramas com **deriva forte** (gradiente de temperatura/solvente), use `lam` menor ($10^5$-$10^6$) e `p` bem pequeno ($10^{-5}$) para que a baseline curve sem cortar os picos.

---

## 🔍 2. Problemas de Detecção (O que é Pico?)

O software está vendo coisas demais ou de menos.

| Sintoma Visual | Diagnóstico | Ação Recomendada | Parâmetro |
| :--- | :--- | :--- | :--- |
| **Muitos picos "fantasmas" (ruído)** | O limite de SNR está muito baixo ou o ruído foi subestimado. | **Aumentar** SNR.<br>**Aumentar** largura mínima. | `snr_threshold` ($\uparrow$)<br>`min_width_seconds` ($\uparrow$) |
| **Picos pequenos reais não detectados** | O software acha que é ruído. | **Diminuir** SNR.<br>**Suavizar** o sinal antes. | `snr_threshold` ($\downarrow$)<br>`sg_window_length` ($\uparrow$) |
| **Topo do pico detectado como 2 picos** | "Ruído serrilhado" no topo do pico (split peak falso). | **Aumentar** a distância mínima ou suavizar. | `min_distance_seconds` ($\uparrow$)<br>`sg_window_length` ($\uparrow$) |
| **Pico do solvente atrapalhando** | O início do cromatograma é uma bagunça. | **Cortar** o início ou filtrar por área. | `t_start_integration`<br>`solvent_rt_cutoff_s` |

---

## 📐 3. Problemas de Integração (Início, Fim e Área)

Os marcadores de início (start) e fim (end) estão nos lugares errados.

| Sintoma Visual | Diagnóstico | Ação Recomendada | Parâmetro |
| :--- | :--- | :--- | :--- |
| **Integração corta a cauda (Tailing)** | O software "acha" que o pico acabou cedo demais (slope alto). | **Diminuir** a sensibilidade do slope. | `slope_threshold_factor` ($\downarrow$)<br>(ex: $0.1 \to 0.05$) |
| **Integração pega muita baseline** | O pico parece "gordo" demais na base. | **Aumentar** a sensibilidade do slope. | `slope_threshold_factor` ($\uparrow$)<br>(ex: $0.1 \to 0.5$) |
| **Áreas inconsistentes (Sobe/Desce)** | A linha de base local está variando muito. | Verificar se a `baseline` global está correta primeiro. | (Volte para a seção 1) |

---

## 🏔️ 4. Problemas de Separação (Picos Colados)

Como o software lida com o vale entre dois picos.

| Sintoma Visual | Diagnóstico | Ação Recomendada | Parâmetro |
| :--- | :--- | :--- | :--- |
| **Usa Drop-line (reto) mas devia ser Skim** | O vale é profundo, mas você quer skimming (pico pequeno no ombro). | **Aumentar** a tolerância de altura do rider. | `height_ratio_rider` ($\uparrow$)<br>(ex: $0.1 \to 0.2$) |
| **Usa Skim mas devia ser Drop-line** | O pico "ombro" é muito grande para ser skimmed. | **Diminuir** a tolerância de altura. | `height_ratio_rider` ($\downarrow$) |
| **Picos fundidos (sem separação)** | O software não viu o vale ou o Rs é alto demais. | **Aumentar** o limiar de Rs para forçar a análise do vale. | `rs_deconv_threshold` ($\uparrow$)<br>(ex: $1.2 \to 1.5$) |
| **Deconvolução falha/ruim** | Picos não gaussianos ou muito ruído. | Tentar forçar Drop-line diminuindo a exigência do vale. | `valley_pct_dropline` ($\downarrow$)<br>(ex: $50 \to 30$) |

---

## 📊 5. Problemas de System Suitability (N, Tf, Rs)

Os números não batem com a validação ou literatura.

| Sintoma Numérico | Causa Provável | Ação |
| :--- | :--- | :--- |
| **N (Pratos) muito baixo** | Suavização excessiva achatou o pico e aumentou a largura ($W_{1/2}$). | **Reduzir** `sg_window_length` (ou zerar para auto). |
| **Tailing Factor (Tf) errado** | O corte da cauda (`end_index`) está errado. | Ajustar `slope_threshold_factor`. Se corta a cauda, Tf melhora artificialmente (erro!). |
| **Resolução (Rs) variando** | A largura na base ($W_{base}$) está instável devido ao ruído. | Usar fórmula EP (meia-altura) é mais robusto que USP (base). |

---

## ⚡ Interações Perigosas ("O Efeito Borboleta")

Mudar um parâmetro pode quebrar outra coisa. Cuidado com estes pares:

1. **Smoothing (`sg_window`) vs. Resolução (`Rs`)**
    * *Ação:* Você aumenta o smoothing para sumir com o ruído.
    * *Efeito Colateral:* Picos próximos se fundem. A altura diminui. A largura aumenta.
    * *Resultado:* **N cai, Rs cai.** Use com moderação.

2. **Baseline (`lam`) vs. Área (`Area`)**
    * *Ação:* Você diminui o `lam` para remover uma "barriga" na baseline.
    * *Efeito Colateral:* A baseline sobe um pouco nas pontas dos picos.
    * *Resultado:* **Área diminui.** A quantificação muda. Mantenha consistente entre padrões e amostras.

3. **Threshold (`snr`) vs. Ruído (`baseline_p`)**
    * *Ação:* Você ajusta o `baseline_p` e a linha de base muda de nível.
    * *Efeito Colateral:* O cálculo de ruído muda. O SNR dos picos muda.
    * *Resultado:* Picos que antes eram detectados podem desaparecer (ou vice-versa) sem você mexer no `snr_threshold`.

---

# 📊 Glossário de Métricas Analíticas (GCAnalyzer)

Este documento detalha todas as métricas calculadas pelo sistema, divididas entre **Métricas Individuais** (características de cada pico), **Métricas Globais** (qualidade da corrida inteira) e **Métricas Comparativas** (alinhamento e similaridade entre corridas).

---

## 1. Métricas Individuais (Por Pico)

Estas informações estão presentes no DataFrame `results_df` e na tabela principal do relatório PDF.

### A. Identificação e Quantificação

| Variável | Nome Completo | Unidade | Descrição e Interpretação |
| :--- | :--- | :--- | :--- |
| **`rt`** | **Tempo de Retenção** | $s$ ou $min$ | Tempo decorrido desde a injeção até o ápice do pico. É a "impressão digital" da identidade do composto. |
| **`rrt`** | **Tempo de Retenção Relativo** | Adimensional | Razão $RT_{pico} / RT_{IS}$ (Padrão Interno). <br>• **Uso:** Corrige flutuações instrumentais (fluxo, temperatura). Permite alinhar corridas diferentes. |
| **`area`** | **Área Integrada** | $u.a. \cdot s$ | Soma da intensidade do sinal sob a curva. Proporcional à quantidade de massa do analito.<br>• Calculada via *trapezoidal* ou *EMG* dependendo do método. |
| **`area_pct`** | **Área Percentual** | $\%$ | Participação do pico na área total integrada da corrida.<br>• **Interpretação:** Pureza relativa (se assumirmos fator de resposta 1:1). |
| **`height`** | **Altura** | $u.a.$ | Intensidade máxima do pico acima da linha de base. |
| **`snr`** | **Signal-to-Noise Ratio** | Adimensional | Razão $Altura / Ruído_{local}$.<br>• **$> 10$:** Quantificável (LOQ).<br>• **$> 3$:** Detectável (LOD).<br>• **$< 3$:** Ruído/Incerto. |

### B. Eficiência e Forma (System Suitability)

Métricas críticas para validar a qualidade da coluna e do método cromatográfico.

| Variável | Nome Completo | Fórmula / Ref | Descrição e Interpretação |
| :--- | :--- | :--- | :--- |
| **`N_plates_ep`** | **Pratos Teóricos (EP)** | $5.54 (tR/W_{1/2})^2$<br>*(Ph. Eur.)* | Mede a eficiência da coluna usando a largura a meia-altura.<br>• **Ideal:** Quanto maior, melhor. Colunas capilares > 10.000.<br>• **Nota:** Mais robusto para picos com cauda que a fórmula USP. |
| **`N_plates_usp`** | **Pratos Teóricos (USP)** | $16 (tR/W_{base})^2$<br>*(USP <621>)* | Fórmula clássica usando a largura na base (tangentes).<br>• Sensível a tailing excessivo. |
| **`tailing_factor_usp`** | **Fator de Cauda (Tf)** | $W_{0.05} / 2f$<br>*(USP)* | Mede a simetria a 5% da altura.<br>• **1.0:** Perfeito.<br>• **$> 1.5$:** Tailing (cauda).<br>• **$< 0.8$:** Fronting (frente).<br>• **Limite aceitável:** Geralmente $0.8 - 2.0$. |
| **`asymmetry_factor_ep`** | **Fator de Assimetria (As)** | $B/A$<br>*(EP)* | Mede a assimetria a 10% da altura.<br>• Similar ao Tf, mas calculado mais acima no pico. |
| **`width_half_s`** | **Largura Meia-Altura** | $s$ | Largura do pico a 50% da altura ($W_{1/2}$). Fundamental para calcular resolução EP. |
| **`sqs`** | **Shape Quality Score** | $0.0 - 1.0$ | Métrica customizada do algoritmo.<br>• **1.0:** Formato Gaussiano perfeito.<br>• **0.0:** Pico distorcido ou ruído. Baseado na razão $\tau/\sigma$ do modelo EMG. |

### C. Resolução e Termodinâmica

Relação do pico com seus vizinhos e com a fase estacionária.

| Variável | Nome Completo | Descrição e Interpretação |
| :--- | :--- | :--- |
| **`Rs_usp`** | **Resolução (USP)** | Separação entre este pico e o **anterior**.<br>• **$\ge 1.5$:** Separação na linha de base (Baseline resolved).<br>• **$< 1.5$:** Coeluição parcial. |
| **`k_prime` ($k'$)** | **Fator de Capacidade** | $(tR - t0) / t0$. Mede o quanto o analito interage com a fase estacionária.<br>• Requer `dead_time_s` configurado.<br>• **Ideal:** $1 < k' < 20$. |
| **`alpha` ($\alpha$)** | **Seletividade** | $k'_2 / k'_1$. Capacidade do sistema químico de distinguir dois compostos.<br>• **1.0:** Coeluição perfeita (sem separação). |

### D. Scoring Composto (CQI)

*Score sintético para facilitar a triagem rápida.*

* **`CQI` (Chromatographic Quality Index):** Média geométrica ponderada (0 a 1) combinando Eficiência ($N$), Resolução ($Rs$), Simetria ($Tf$) e Ruído ($SNR$).
  * **Uso:** Ordenar picos do "melhor" para o "pior" para auditoria rápida.

---

## 2. Métricas Globais (Da Corrida)

Estas métricas avaliam a saúde da análise como um todo (instrumento + amostra).

| Variável | Descrição e Interpretação |
| :--- | :--- |
| **`baseline_drift`** | Variação absoluta da linha de base (Início vs. Fim).<br>• **Alto:** Indica gradiente de temperatura agressivo ou sangramento da coluna (bleeding). |
| **`global_snr`** | Razão entre o pico mais alto da corrida e o ruído de fundo global.<br>• Indica a sensibilidade geral da injeção. |
| **`total_integrated_area`** | Soma das áreas de todos os picos.<br>• **Uso:** Se variar muito entre injeções do mesmo vial, indica problema no injetor ou vazamento. |
| **`area_coverage_pct`** | Qual porcentagem do sinal total (acima da baseline) foi reconhecida como pico.<br>• **Baixo (< 80%):** Muito ruído não integrado ou "morro" de linha de base não resolvido. |
| **`OQS` (Overall Quality Score)** | Nota geral da corrida (0 a 1). Combina `CQI` médio, estabilidade da baseline e SNR global. |

---

## 3. Métricas de Comparação (Fingerprinting & Alinhamento)

Usadas quando se compara duas corridas (ex: `Amostra vs Padrão` ou `Lote A vs Lote B`).

### A. Alinhamento de Picos (Binning)

O sistema tenta casar picos entre corridas.

* **`rt_shift` (Deslocamento de RT):** Diferença temporal ($RT_{amostra} - RT_{ref}$) para o mesmo pico.
  * **Causa:** Variação de fluxo, temperatura ou envelhecimento da coluna.
  * **Correção:** O uso de **RRT** (via Padrão Interno) minimiza o impacto visual desse shift.
* **`rrt_cv_pct`:** Coeficiente de Variação do RRT para um mesmo pico em múltiplas corridas.
  * **Ideal:** $< 0.5\%$. Se alto, o alinhamento está falhando.

### B. Comparação de Sinal Total (Fingerprint)

Compara o cromatograma inteiro ponto a ponto, não apenas a lista de picos.

| Variável | Descrição | Interpretação |
| :--- | :--- | :--- |
| **`pearson_r`** | Correlação Linear (0 a 1). | **> 0.990:** Perfis praticamente idênticos.<br>**< 0.900:** Diferenças significativas (impurezas, degradação). |
| **`cosine_similarity`** | Similaridade de Cosseno. | Similar ao Pearson, mas insensível à magnitude (intensidade). Foca na "forma" do perfil. |
| **`spectral_contrast_angle`** | Ângulo entre vetores. | **0°:** Idêntico. **90°:** Totalmente diferente. Usado em bibliotecas espectrais. |
| **`verdict`** | Julgamento automático. | `IDENTICAL`, `SIMILAR`, `ACCEPTABLE`, `DIFFERENT`. Baseado nos limiares de Pearson. |

---

## 4. Tipos de Integração (`integration_method`)

Explica *como* a área foi calculada. Aparece no relatório e nos gráficos.

1. **`EMG`:** Ajuste matemático (Exponentially Modified Gaussian). O mais preciso para picos assimétricos. Separa o sinal do pico do ruído/cauda.
2. **`TRAPEZOID`:** Soma simples dos pontos (regra do trapézio). Usado quando o ajuste EMG falha ou diverge. É o método "clássico" robusto.
3. **`DROP_LINE`:** Separação de picos colados por uma linha vertical no vale. Usado quando o vale é profundo (> 50%).
4. **`TANGENT_SKIM`:** Separação de um pico pequeno ("rider") no ombro de um grande ("parent") por uma linha tangente.
5. **`DECONVOLUTION`:** Separação matemática de múltiplos picos sobrepostos (ajuste multi-pico). Usado para picos complexos/misturados.
6. **`FORCED`:** Integração manual forçada pelo usuário via `force_integration_windows`. Ignora filtros de qualidade.

---

# 📐 Modelos Matemáticos e Parâmetros de Forma (Gaussian vs. EMG)

O **GCAnalyzer** não apenas "soma" a área sob o gráfico. Ele tenta ajustar modelos matemáticos aos dados brutos para entender a física por trás de cada pico.

Isso é fundamental para separar picos sobrepostos (deconvolução) e distinguir o sinal real do ruído.

---

## 1. O Modelo Gaussiano (Ideal)

Na teoria cromatográfica perfeita (difusão pura), todo pico é uma **Gaussiana**. É uma curva em forma de sino, perfeitamente simétrica.

### Fórmula

$$f(x) = A \cdot \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)$$

### Parâmetros

| Parâmetro | Nome | O que controla? | Interpretação Cromatográfica |
| :--- | :--- | :--- | :--- |
| **$\mu$ (Mu)** | Média | **Posição** | O tempo de retenção exato do ápice ($t_R$). |
| **$\sigma$ (Sigma)** | Desvio Padrão | **Largura** | A dispersão da banda. Quanto maior o $\sigma$, mais "gordo" e baixo é o pico (menor eficiência $N$). |
| **$A$** | Amplitude | **Altura/Área** | A intensidade máxima (se normalizado, relaciona-se à área). |

### Quando ocorre?

* Colunas novas e muito eficientes.
* Analitos não polares.
* Picos muito rápidos (iniciais).

---

## 2. O Modelo EMG (Realidade)

Na prática, picos cromatográficos quase sempre têm uma "cauda" (tailing). Isso ocorre devido a efeitos de transferência de massa, adsorção na coluna ou volume morto no detector.

Para modelar isso, usamos a **Exponentially Modified Gaussian (EMG)**. É a soma matemática (convolução) de uma Gaussiana (o pico ideal) com um Decaimento Exponencial (o atraso do sistema).

### Parâmetros Extras

Além de $\mu$ e $\sigma$, o EMG introduz:

| Parâmetro | Nome | O que controla? | Interpretação Cromatográfica |
| :--- | :--- | :--- | :--- |
| **$\tau$ (Tau)** | Relaxamento | **Assimetria (Cauda)** | O tempo médio que as moléculas ficam "presas" ou atrasadas no sistema. |
| **$\lambda$** | Taxa | **Inverso de Tau** | $\lambda = 1/\tau$. (Menos usado nos relatórios, mas matematicamente equivalente). |

---

## 3. Entendendo a Relação $\tau$ e $\sigma$

A "personalidade" de um pico é definida pela briga entre a difusão ($\sigma$) e o atraso ($\tau$).

### O Ratio $\tau / \sigma$ (Tau-Sigma Ratio)

Esta é a métrica mais importante para o diagnóstico de forma no código (`shape_quality_score`).

1. **$\tau \approx 0$ (ou $\tau \ll \sigma$):**
    * **Forma:** Gaussiana Pura. Simétrico.
    * **Diagnóstico:** Sistema ideal.

2. **$\tau \approx \sigma$:**
    * **Forma:** Levemente assimétrico. O lado direito desce mais devagar que o esquerdo sobe.
    * **Diagnóstico:** Típico de análises reais aceitáveis.

3. **$\tau \gg \sigma$ (Tau muito maior que Sigma):**
    * **Forma:** "Dente de Tubarão". O pico sobe rápido e desce muito devagar (cauda longa).
    * **Diagnóstico:** Problema severo. Sítios ativos na coluna, volume morto no injetor ou degradação da fase estacionária.

---

## 4. Visualização dos Parâmetros

Imagine alterar apenas um parâmetro enquanto mantém os outros fixos:

### Aumentando apenas $\sigma$ (Sigma)

* O pico fica **mais largo** simetricamente.
* O ápice desce.
* A área se mantém (se $A$ for compensado).
* **Causa:** Difusão longitudinal (fluxo muito baixo) ou caminho muito longo.

### Aumentando apenas $\tau$ (Tau)

* O início do pico (subida) quase não muda.
* O ápice se desloca levemente para a direita.
* O final do pico (descida) se estica muito.
* **Causa:** Adsorção química ("column drag").

---

## 5. Métricas Derivadas no `results_df`

O GCAnalyzer usa $\sigma$ e $\tau$ para calcular métricas de qualidade (SQS) e converter para normas farmacopeicas.

### A. Shape Quality Score (SQS)

Calculado como:
$$SQS = \exp\left(-\max\left(\frac{\tau}{\sigma} - 1, 0\right)\right)$$

* **Se $\tau < \sigma$:** SQS = 1.0 (Qualidade Perfeita/Gaussiana).
* **Se $\tau = 2\sigma$:** SQS cai rapidamente.
* **Interpretação:** Um SQS baixo indica que o pico não é confiável para cálculos de eficiência padrão (N) baseados em modelos Gaussianos.

### B. Tailing Factor (Tf) vs. $\tau$

Embora o USP Tailing Factor ($T_f$) seja medido geometricamente ($W_{0.05}/2f$), ele é fortemente correlacionado com $\tau$:

* **Gaussiana ($\tau=0$):** $T_f = 1.0$
* **EMG ($\tau > 0$):** $T_f > 1.0$

> **Nota:** Se o ajuste EMG falhar (não convergir), o parâmetro `A_param` será `NaN` e o software reverterá para integração trapezoidal (área geométrica simples). Nesses casos, $\sigma$ e $\tau$ não estarão disponíveis.

---

## 6. Resumo Matemático para Debugging

Se você vir estes valores nos logs ou no CSV exportado:

| Cenário | Valores Observados | O que significa? |
| :--- | :--- | :--- |
| **Pico Perfeito** | $\sigma > 0$, $\tau \approx 0$, $A > 0$ | Ajuste Gaussiano bem sucedido. |
| **Pico Real** | $\sigma > 0$, $\tau > 0$, $A > 0$ | Ajuste EMG bem sucedido. |
| **Falha de Ajuste** | $\sigma = \text{NaN}$, $\tau = \text{NaN}$ | O algoritmo `curve_fit` não conseguiu ajustar a curva. Área calculada por trapézio. |
| **Ruído?** | $\sigma$ muito pequeno (< 0.05s) | Provavelmente um "spike" eletrônico, não um pico cromatográfico real. |
| **Deconvolução** | Múltiplos sets de ($\mu, \sigma, \tau$) | O algoritmo separou matematicamente picos que estavam somados. |

---

# 📐 Fórmulas Cromatográficas: USP vs. EP

Este documento detalha as equações matemáticas utilizadas pelo `GCAnalyzer` para calcular eficiência, resolução e simetria.

O software implementa **ambas** as normas (USP e EP) simultaneamente. É crucial entender a diferença, pois elas medem o pico em alturas diferentes, gerando resultados numéricos distintos para o mesmo pico.

---

## 1. Eficiência da Coluna (N - Pratos Teóricos)

Mede o quão "fino" é o pico em relação ao seu tempo de retenção. Quanto maior o $N$, menor o alargamento da banda e melhor a performance da coluna.

### A. Fórmula USP (United States Pharmacopeia)

Historicamente baseada na largura da base (tangentes).

$$ N_{USP} = 16 \left( \frac{t_R}{W_{base}} \right)^2 $$

* **$t_R$**: Tempo de retenção do ápice.
* **$W_{base}$**: Largura do pico na base (estimada a **5% da altura** no software para robustez digital).
* **Características:** Muito sensível ao alargamento na base (tailing). Se o pico tiver cauda, o $W_{base}$ aumenta muito e o $N_{USP}$ cai drasticamente.

### B. Fórmula EP (European Pharmacopoeia)

Baseada na largura a meia-altura (FWHM).

$$ N_{EP} = 5.54 \left( \frac{t_R}{W_{1/2}} \right)^2 $$

* **$W_{1/2}$**: Largura do pico a **50% da altura**.
* **Características:** Mais robusta e reproduzível para picos assimétricos, pois ignora a base larga. Geralmente resulta em valores de $N$ ligeiramente maiores que a USP em picos reais.

---

## 2. Simetria do Pico (Tailing vs. Asymmetry)

Aqui reside a maior confusão em cromatografia. USP e EP usam nomes diferentes e alturas de medição diferentes.

### A. USP Tailing Factor ($T_f$)

Métrica padrão exigida pela FDA/USP. Medida a **5% da altura**.

$$ T_f = \frac{W_{0.05}}{2 \cdot f} $$

* **$W_{0.05}$**: Largura total do pico a 5% da altura.
* **$f$**: Distância da borda frontal (subida) até a projeção do ápice (tempo de retenção), medida a 5% da altura.
* **Interpretação:**
  * $T_f = 1.0$: Perfeitamente simétrico.
  * $T_f > 1.0$: Tailing (Cauda). $T_f \le 2.0$ é o limite comum.
  * $T_f < 1.0$: Fronting (Frente).

### B. EP Asymmetry Factor ($A_s$)

Métrica padrão na Europa. Medida a **10% da altura**.

$$ A_s = \frac{b_{0.1}}{a_{0.1}} $$

* **$b_{0.1}$**: Distância da projeção do ápice até a borda traseira (cauda), a 10% da altura.
* **$a_{0.1}$**: Distância da borda frontal (subida) até a projeção do ápice, a 10% da altura.
* **Diferença:** Como é medido mais alto no pico (10% vs 5%), o $A_s$ tende a ser numericamente menor que o $T_f$ para o mesmo pico com cauda.

> **Nota do Software:** No relatório, `tailing_factor_usp` usa a regra dos 5% e `asymmetry_factor_ep` usa a regra dos 10%.

---

## 3. Resolução ($R_s$)

Mede a separação entre dois picos adjacentes (Pico 1 e Pico 2). O cálculo é feito sempre em relação ao pico **anterior**.

### A. Fórmula USP (Tangente/Base)

$$ R_{s, USP} = \frac{2 (t_{R2} - t_{R1})}{W_{base1} + W_{base2}} $$

* Usa as larguras na base ($W_{base}$).
* **Critério:** $R_s \ge 1.5$ indica separação na linha de base (erro < 1%).
* **Sensibilidade:** Penaliza fortemente picos com cauda, pois $W_{base}$ inclui a cauda. É a métrica mais conservadora ("pior caso").

### B. Fórmula EP (Meia-Altura)

$$ R_{s, EP} = \frac{1.18 (t_{R2} - t_{R1})}{W_{1/2, 1} + W_{1/2, 2}} $$

* Usa as larguras a meia-altura ($W_{1/2}$).
* **Diferença:** Se os picos tiverem muita cauda, o $R_{s, EP}$ dará um valor maior (mais otimista) que o USP, pois ignora a sobreposição que ocorre lá embaixo na base.

---

## 4. Parâmetros Termodinâmicos

Estes dependem da configuração correta do Tempo Morto ($t_0$ ou `dead_time_s`). Se $t_0 = 0$, estes valores não são calculados.

### A. Fator de Retenção / Capacidade ($k'$)

Mede quanto tempo o analito fica "preso" na fase estacionária em relação à fase móvel.

$$ k' = \frac{t_R - t_0}{t_0} $$

* **$t_0$**: Tempo morto (tempo para um composto não retido atravessar a coluna).
* **Interpretação:**
  * $k' = 0$: O composto não interagiu (saiu no tempo morto).
  * $k' < 2$: Saiu muito rápido, risco de interferência da matriz.
  * $k' > 20$: Tempo de análise excessivamente longo.

### B. Seletividade ($\alpha$)

Mede a capacidade química da coluna de distinguir dois compostos (Pico 2 vs Pico 1).

$$ \alpha = \frac{k'_2}{k'_1} = \frac{t_{R2} - t_0}{t_{R1} - t_0} $$

* **Interpretação:**
  * $\alpha = 1.0$: Coeluição perfeita (impossível separar, não importa a eficiência $N$).
  * $\alpha > 1.1$: Separação fácil.

---

## 5. Signal-to-Noise Ratio (SNR)

O `GCAnalyzer` utiliza uma abordagem estatística robusta para o ruído, conforme diretrizes da ASTM/USP modernas.

$$ SNR = \frac{Height}{Noise_{sigma}} $$

* **Height:** Altura do pico acima da linha de base interpolada localmente.
* **Noise ($ \sigma $):** Calculado pelo método **MAD (Median Absolute Deviation)** na região do pico ou globalmente.
  * Fórmula do ruído: $\sigma \approx 1.4826 \times \text{median}(|x - \text{median}(x)|)$
  * Isso evita que "spikes" isolados inflem o valor do ruído, fornecendo um SNR mais realista.

---

## Resumo Visual de Medição

Para garantir que você está olhando a métrica certa:

| Altura de Medição | Métrica USP | Métrica EP | Símbolo no Código |
| :--- | :--- | :--- | :--- |
| **50% (Meia)** | — | Eficiência ($N$) | `W_half_s` |
| **10%** | — | Assimetria ($A_s$) | `asymmetry_factor_ep` |
| **5% (Base)** | Eficiência ($N$), Tailing ($T_f$) | — | `W_base_s` |

> **Atenção:** Embora a USP tradicionalmente use tangentes, métodos computacionais modernos (como este software) aproximam a "base" pela largura a 5% ou 10% dependendo da configuração, pois tangentes são matematicamente instáveis em sinais com ruído digital.
