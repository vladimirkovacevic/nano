
# Important for all figures:
Create figures and their captions for the advance materials journal.  Start with notebook @notebooks/nano_figures.ipynb and modify it (extend it). Reuse existing code for visualizing vulcano plots, gene heatmaps and distribution of enriched pathways. Sample 1 is treated with 40nm PSNPs, so we call it 40nm, Sample 2 - 200nm PSNPs and Sample 3 with combination 40+200nm, which is why we call it 40+200nm. Last Sample 4 is called Untreated. In the notebook add figures in the given order and according to a given description. Define color palete at the top and keep it consistent across all figures. Save all figures in fig folder. At the end of each figure generate latex for placing all figures together with their caption and notify me to review it.
Each figure should be consisted of subplots (subfigures) a, b, c, d... They should be one figure, don't separate a and b into one and c and d into the other figure.
Don't put text as figure, but try to think of a way to construct this figure with that given text.
Don't close figures let user see them in notebooks/nano_figures.ipynb.
All text must be large enough so it is readable (as large as possible). Use Arial font everywhere!
Generate only PDF of the figures, not PNGs.
Execute, test and debug notebook!
Professional code with minimal complexity, and a lot of comments, for generating main and supplementary figures save in notebooks/nano_figures.ipynb


DEGs: results/de_results_CoDi_dist.xlsx
GSE: data/Pathways_KEGG_cleaned.xlsx (Sheet names contain Cell type and condition (sample), just remove spaces when decoding them)
Data: data/merged.h5ad (use CoDi_dist cell type annotation)


# Figure 1 - Overview
## Figure 1A Microfluidic chip design, blood flow, PSNP exposure (CTRL, 40 nm, 200 nm, 40 + 200 nm; exposure 24 h) → scRNA-seq. Explicitly label “24 h exposure” and “flow” on the scheme (even with arrows)
## Fig ​1B,C - UMAP of all immune cells across conditions
(Data-driven, global view)
What it should show
UMAP embedding of all cells
Make 1B = colored by cell type
Make 1C = colored by condition
Why it matters
Demonstrates:
Successful scRNA-seq
Clear immune cell identities
No batch collapse or technical artifacts
Fig. 1D - Immune cell composition across PSNP conditions
(Quantitative validation)
What it should show
Bar plot or stacked bar plot
X-axis: Experimental condition (Control, 40 nm, 200 nm, 40 + 200 nm)
Y-axis: Percentage of total cells
Colors: cell types (same colors as Fig. 1B)
Why it matters
Shows:
Cell populations are largely preserved
DEG results are not driven by compositional shift


# Figure 2. Size-dependent transcriptional impact of PSNP exposure across immune cell types
## Fig. 2A  DEG burden per cell type and experimental condition
(Primary quantitative panel): Shows bar plot or grouped bar plot
X-axis: Immune cell types
 (B cells, CD4⁺ T, Cytotoxic T, Monocytes, NK)
Y-axis: Number of DEGs (adjusted p < 0.05)
Color-coded bars:
40 nm PSNPs
200 nm PSNPs
40 + 200 nm PSNPs
Expected: 40 nm → strongest transcriptional impact; our data 200 nm more DEGs
Mixed exposure ≠ simple sum
Design tips (Advanced Materials)
Use side-by-side bars per cell type
Add small numeric labels above bars (optional but powerful)
Results text it supports: 200 nm → more DEGs. This is actually interesting and non-intuitive → good
## Figure 2B DEG burden normalized to cell abundance 
(calculate DEG density per 1,000 cells)
use data/merged.h5ad (use CoDi_dist cell type annotation)
Step-by-step:
Step 1 — Get cell counts
For each cell type × condition, count the number of cells used in the DEG test.

Example:
Cell type
Condition
- cells
Monocytes
40 nm
2,400
B cells
40 nm
4,800
NK
40 nm
350

Step 2 — Get DEG counts
From your DEG table (adjusted p < 0.05, same threshold everywhere):
Cell type
Condition
- DEGs
Monocytes
40 nm
1,200
B cells
40 nm
600
NK
40 nm
150


Step 3 — Normalize
Compute:
DEG density=Number of DEGsNumber of cells×1000\textbf{DEG density} = \frac{\text{Number of DEGs}}{\text{Number of cells}} \times 1000DEG density=Number of cellsNumber of DEGs​×1000

Example:
Cell type
Condition
DEG density
Monocytes
40 nm
500
B cells
40 nm
125
NK
40 nm
429


Step 4 — Plot it
Plot DEG density, not raw DEG counts.
Plot design:
X-axis: Cell type
Y-axis: DEGs per 1,000 cells
Color: PSNP condition (40 nm, 200 nm, mix)
Error bars: Optional (bootstrapped CI if you want to be fancy)

How to interpret (example text)
“When normalized to cell abundance, monocytes and NK cells displayed a disproportionately high DEG density, indicating enhanced transcriptional sensitivity to PSNP exposure.”


## Figure 2C Directionality of transcriptional changes

Stacked bar plot
X-axis: Cell types
Y-axis: Percentage of DEGs
Stacked segments:
Upregulated genes
Downregulated genes
Separate stacks for each condition (faceted or grouped)
We are observing early adaptation, which can be seing by maintain or upregulate core biosynthetic machinery while globally repressing non-essential transcriptional programs, in light with ribosome/translation signal, the absence of strong inflammatory pathways, the short exposure (24 h), the microfluidic (physiological) context.

## Figure 2D  Full volcano pots for all cell types × all conditions
No thresholds on x-axis
No gene labels


# FIGURE 3. Conserved transcriptional programs induced by PSNP exposure

Figure 3 establishes the first biological principle of the paper: PSNP exposure induces a conserved transcriptional program across immune cell types, dominated by translation and RNA related processes.
Increase maximally x-ticks and y-ticks

## Figure 3A: KEGG pathways
Pathway × Cell type heatmap (Arial - size large - big)
Rows: selected KEGG pathways (below) with NES scores taken from data/Pathways_KEGG_cleaned.xlsx
Columns: Cell type × condition (Monocite 40, Monocite 200, Monocite Mix, B-cells 40nm, B-cells 200, …)
Color: NES (red ↑, blue ↓), Grey: no enrichment
Increase size of x-ticks and y-ticks and show x-ticks for all samples and all conditions
Show on y-axis Pathways (start with "-") gathered into groups (color square):
🟦 Translation & RNA regulation 
- Ribosome
🟩 Mitochondrial metabolism & bioenergetics
- Oxidative phosphorylation
- Citrate cycle (TCA cycle)

- Pentose phosphate pathway
🟥 Innate immune sensing & inflammatory signaling
- TNF signaling pathway
- IL-17 signaling pathway
🟪 Antigen presentation & immune communication
- Antigen processing and presentation
(color font of each pathway with the color of the group and put groups in the legend)
## Figure 3B:
Split into three clearly visible blocks:
Block 1 — Translation & RNA regulation 
RPL10
RPL11
RPL14
RPL24
RPS2
RPS12
RPS27A
RPS28
RPS3
Block 2 — Antigen processing & immune communication
Genes:
HLA-B,
HLA-C,
HLA-DMA,
HLA-DPA1,
CD74.
Block 3 — Minimal lineage markers (context, not mechanisms - not for interpretation, only to help read the heatmap). 
Genes:
MS4A1, CD79A
CD247, GRAP2

Structure
Rows: genes (grouped by block, not scattered across blocks)
Columns: Cell type × condition
Values: log2FC

Plot type: dot plot (Seurat-style or equivalent)
Y-axis: Blocks, Genes
X-axis: Cell type × condition
Group by cell type
Within each group: Control | 40 nm | 200 nm | 40+200 nm
Encoding
Dot color: average log2 fold change (vs control)
Make font size as large enough to it can be read from any computer
--------------------------------------------------------

# FIGURE 4  Size-dependent transcriptional remodeling in monocytes following PSNP exposure

## Fig. 4A Size-dependent differential gene expression in monocytes
Three volcano plots (side by side) with all visible outliers - no cut off!
Left: CD14⁺ monocytes, 40 nm vs control
Midle: CD14⁺ monocytes, 200 nm vs control
Right: CD14⁺ monocytes, 40+200 nm vs control


Axes (identical for both plots)
X-axis: log2 fold change
Y-axis: −log10(adjusted p-value)
Genes highlighted (verified DEGs)
One volcano per condition, label only anchor genes (small black circle, not star as marker and gene names with light gray background. Organize labels so they don't overlap between each other):
40 nm: CD36, FGL2, PDLIM4, GPX3 
200 nm: SPP1, GPC4, AREG, ENHO 
40+200: GPC4, FGL2, MRC1L1, SIGLEC7 


## Fig. 4B  Pathway-anchored gene heatmap. Purpose: show state transitions, not overlap.
Rows (y-axis) = should show both module and anchor genes separated by comma (grouped by module in provided order)
Columns = conditions (they can be narrow, to leave enough space for y-ticks label containing both module name and gene name)
Y-axis modules: Transcriptional / RNA control, Vesicle trafficking, Metabolism Regulatory signaling

Here are modules and their genes:
🟦 MODULE A —Translation & RNA regulation (40 nm–specific)
Pathways: Ribosome, Spliceosome, RNA transport
Anchor genes: RPL41, RPL7, RPS16, FAU

🟦 MODULE B — Vesicle trafficking & intracellular transport (200 nm–specific)
Pathways: Synaptic vesicle cycle
Anchor genes: GPC4, SPP1, AREG 

🟦 MODULE C — Mitochondrial metabolism & bioenergetics (200 nm + 40+200)
Pathways: Fatty acid elongation, Citrate cycle, Pentose phosphate pathway
Anchor genes: CD36, ENHO, GPX3 

🟦 MODULE D — Innate immune sensing and controlled signaling(strongest in 40+200)
Pathways: TNF signaling, IL-17 signaling, Viral protein–cytokine interaction
Anchor genes: FGL2, MRC1L1, SIGLEC7 

## Figure 4C Coherent pathway modulation in monocytes is size dependent
Dot plot (monocytes only)
X-axis: PSNP condition (40 nm, 200 nm, 40+200 nm)
Y-axis: Module, Pathway (Important: Add all pathways from all modules below!)
Modules (pathways are given in brackets):
- Translation & RNA regulation (Ribosome, Spliceosome, mRNA surveillance, RNA transport)
- Vesicle trafficking & intracellular transport (Synaptic vesicle cycle)
- Mitochondrial metabolism & bioenergetics (Citrate cycle (TCA), Pentose phosphate pathway, Fatty acid elongation, Thermogenesis)
- Innate immune sensing and controlled signaling(TNF signaling pathway, IL-17 signaling pathway, Viral protein–cytokine receptor interaction)
(Assign color to each of 4 modules and use it for font of its pathways. Put colors with names of modules in legend)



## Figure 4D
Blank - add only placeholder with sign - GOOGLE DRAWING
change wording based on newly chosen genes:
Left panel (40 nm): “Transcriptional tuning via translation and RNA regulation ”
Middle panel (200 nm): “Membrane dynamics & extracellular interface remodeling”
Right panel (40+200): “Integrated metabolic–regulatory innate immune state (non-additive)”


# FIGURE 5 — Adaptive immune cells display restrained and lineage-specific transcriptional responses to PSNP exposure
## FIGURE 5A — Volcano plots per condition

label ONLY these genes (small black circle, not star as marker and gene names with light gray background. Avoid overlaping of gene labels):

B cells genes:
40 nm:  MS4A7, LTA

200 nm: FCER2 (CD23), DOK3

40 + 200 nm: FCER2, TNF, HLA-C


CD4⁺ T cells
40 nm:  CCR7, IL7R 
200 nm: CCR7, LTB

40 + 200 nm: CCR7, IFITM1 


## FIGURE 5B — Pathway-anchored heatmap (adaptive immunity)
Show name of the module and gene name on y-axes:
🟦 MODULE A — Translation & RNA regulation (here is minimal in contrast to Figure 4, only tuning), so it can be separate with a line, or in liter color)
RPL7, RSP16, FAU

 
 
🟦 MODULE B — Antigen processing & immune communication 
B cells
TNF
DOK3
HLA-C
CD4⁺ T cells
IFITM1


🟦 MODULE C —  Adaptive immune receptor signaling 

B cells
FCER2
MS4A7
T cells
CCR7
IL7R



🟦 MODULE C —  Cell adhesion, migration & tissue interaction
CCR7
ICAM2
FCER2

 # FIGURE 6 — Innate vs adaptive immune divergence under PSNP exposure
## Figure 6A — Cross-cell-type gene response summary (data-driven)
Purpose: to compare innate and adaptive immune cells directly, using the same representative genes, without re-showing volcanoes or heatmaps.
This panel bridges:
Fig. 4 (monocytes)
Fig. 5 (B/CD4)
What to show: A dot plot (Seurat-style), very clean.
X-axis:
Cell types:
Monocytes
B cells
CD4⁺ T cells
Y-axis:
Representative genes (already validated in earlier figures)
Gene selection 
🟦 Innate (monocyte) anchor genes
(from Fig. 4, red color used in other figures used for font size of genes)
RPL7
RPL41 
CD36 
ENHO 
FGL2
SIGLEC7 –



🟩 Adaptive (B/CD4) anchor genes
(from Fig. 5, , DEG ∩ KEGG verified, blue color used in other figure used for font size of genes)
FCER2 
CCR7 
IL7R 
HLA-C 
ICAM2 

Create legend with Inate and Adaptive colors

Encoding
Dot color: avg log2FC (vs control)
Dot size: % of cells expressing gene (calculated from merged.h5ad)
What the reader should see instantly: 
Strong, condition-specific modulation in monocytes
Weaker, parallel modulation in adaptive cells
Minimal cross-over between gene sets
This visually proves innate vs adaptive divergence.

## Figure 6B — State coherence vs transcriptional modulation
Purpose: to elevate the discussion from genes to response logic.
What to show: a schematic comparison, two side-by-side blocks:
Left: Innate immunity (Monocytes)
Label:
“Pathway-coherent state reprogramming”
Bullet icons:
Translation & RNA regulation (40 nm)
Vesicle trafficking (200 nm)
Metabolic-regulatory signaling  (40+200)
Arrow indicating non-additive behavior
Right: Adaptive immunity (B + CD4⁺)
Label:
“Distributed transcriptional tuning”
Bullet icons:
Immune communication 
Receptor regulation 
Cell positioning
Explicit text:
“No state transition”
“Limited non-additivity”

## Figure 6C — Conceptual model: immune system response to PSNP exposure
Purpose: this is the take-home visual of the entire manuscript.
What to show (BioRender-style)
A horizontal flow:
PSNP exposure (size & complexity)
          ↓
Innate immune sensing (monocytes)
  → size-dependent state reprogramming
  → non-additive responses
          ↓
Adaptive immune modulation
  → transcriptional tuning
  → preserved identity

Visual elements
PSNPs (40 nm, 200 nm, mixed)
Monocyte icon (central, emphasized)
B cell / CD4 icon (secondary)
Arrows showing information flow, not activation


Text inside the panel (minimal)

“Physiological exposure context (microfluidics)”
“Physiological exposure favors controlled immune adaptation over inflammatory activation.”


# Supplementary Figure 1 scRNA-seq quality control and robustness
Purpose:To demonstrate that all downstream DEG and pathway analyses are technically sound and not driven by artifacts. Answering on question: “How confident are you in the scRNA-seq data quality?”
This figure is consisted of 4 panels (A, B, C, D): S1A, S1B, S1C and S1D


## Supplementary Figure S1A  Per-cell QC metrics
Use data from merged.h5ad to measure QC according to best practices for scRNA. Use CoDi_dist annotation from adata.obs
What to show
Violin or box plots of:
Number of detected genes per cell
Total UMIs per cell
% mitochondrial reads
Grouping
Facet by condition
Color by cell type (optional)
Why: shows:
No RNA degradation
No systematic shifts across conditions
Supports interpretation of dominant downregulation
This is your defense against Fig. 2C criticism.
## Supplementary Figure S1B QC metrics per cell type
What to show
Same metrics as 1A
Aggregated per cell type
Why: shows that:
Rare populations (NK, CD8⁺) are not low-quality
Differences in DEG burden are biological, not technical
## Supplementary Figure S1C — Pseudobulk PCA (QC, not biology;  Do not interpret this in the Results.
 Mention it only briefly as QC)
What to show
PCA of pseudobulk transcriptomes
One point = one sample (cell type × condition)
Color
PSNP condition
Why: shows no batch effects, modest separation by particle size (if present)
## Supplementary Figure S1D Marker gene validation
What to show
Dot plot for canonical markers per cell type:
Monocytes: LST1, S100A8/9
B cells: MS4A1, CD79A
T cells: CD3D/E
NK: NKG7
Why: confirms correct annotation; prevents reviewer doubts


# SUPPLEMENTARY FIGURE S2  Analytical support for Figures 2 and 3
Supplementary Fig. S2 provides completeness and robustness for: DEG burden (Fig. 2), Directionality (Fig. 2C), Program identification (Fig. 3)


## Fig. S2A Full volcano plots (supports: Fig. 2A, Fig. 2C, Fig. 3C)
What to show
Volcano plots for all cell types × all conditions
Same thresholds
Same axis limits
Minimal or no gene labels
Limit x-axis [-6,6]

## Fig. S2B Full GO / KEGG enrichment results (supports Fig. 3 pathway interpretation)
Dot plots or bar plots of all enriched GO / KEGG terms
One panel per cell type



# SUPPLEMENTARY FIGURE S3 Full GO / KEGG analysis What to show
Dot plots of all enriched KEGG terms, one panel per cell type.
Use data/Pathways_KEGG_cleaned.xlsx
Use NES values.
y-axis: Pathways
x-axis: Condition (40nm, 200nm, 40+200nm)
.


# SUPPLEMENTARY FIGURE S4
UMAP for MRPL41 gene
4 Panels: Control | 40nm | 200nm | 40+200nm in Monocytes
