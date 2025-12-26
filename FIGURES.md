Important:
Create figures and their captions for the advance materials journal.  Start with notebook @notebooks/nano_figures.ipynb and modify it (extend it). Reuse existing code for visualizing vulcano plots, gene heatmaps and distribution of enriched pathways. Sample 1 is treated with 40nm PSNPs, so we call it 40nm, Sample 2 - 200nm PSNPs and Sample 3 with combination 40+200nm, which is why we call it 40+200nm. Last Sample 4 is called Untreated. In the notebook add figures in the given order and according to a given description. Define color palete at the top and keep it consistent across all figures. Save all figures in fig folder. At the end of each figure generate latex for placing all figures together with their caption and notify me to review it.
Each figure should be consisted of subplots (subfigures) a, b, c, d... They should be one figure, don't separate a and b into one and c and d into the other figure.
Don't put text as figure, but try to think of a way to construct this figure with that given text.
Don't close figures let user see them in jupyter notebook.
Execute, test and debug notebook!



Figure 1 – Platform and experimental design
Fig. 1A Microfluidic chip design, blood flow, PSNP exposure (CTRL, 40 nm, 200 nm, 40 + 200 nm; exposure 24 h) → scRNA-seq. Explicitly label “24 h exposure” and “flow” on the scheme (even with arrows)
Fig ​1B,C - UMAP of all immune cells across conditions
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
Supplementary Fig. 1 scRNA-seq quality control and robustness
Purpose:To demonstrate that all downstream DEG and pathway analyses are technically sound and not driven by artifacts. Answering on question: “How confident are you in the scRNA-seq data quality?”
Fig. S1A  Per-cell QC metrics
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
Fig. S1B QC metrics per cell type
What to show
Same metrics as 1A
Aggregated per cell type
Why: shows that:
Rare populations (NK, CD8⁺) are not low-quality
Differences in DEG burden are biological, not technical
Fig. S1C — Pseudobulk PCA (QC, not biology;  Do not interpret this in the Results.
 Mention it only briefly as QC)
What to show
PCA of pseudobulk transcriptomes
One point = one sample (cell type × condition)
Color
PSNP condition
Why: shows no batch effects, modest separation by particle size (if present)
Fig. S1D Marker gene validation
What to show
Dot plot for canonical markers per cell type:
Monocytes: LST1, S100A8/9
B cells: MS4A1, CD79A
T cells: CD3D/E
NK: NKG7
Why: confirms correct annotation; prevents reviewer doubts

Figure 2. Size-dependent transcriptional impact of PSNP exposure across immune cell types
Fig. 2A  DEG burden per cell type and experimental condition
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
Fig. 2B DEG burden normalized to cell abundance 
(calculate DEG density per 1,000 cells)
Step-by-step
Step 1 — Get cell counts
For each cell type × condition, count the number of cells used in the DEG test.Example:
Cell type
Condition
# cells
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
# DEGs
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

Fig. 2C Directionality of transcriptional changes

Stacked bar plot
X-axis: Cell types
Y-axis: Percentage of DEGs
Stacked segments:
Upregulated genes
Downregulated genes
Separate stacks for each condition (faceted or grouped)
We are observing early adaptation, which can be seing by maintain or upregulate core biosynthetic machinery while globally repressing non-essential transcriptional programs, in light with ribosome/translation signal, the absence of strong inflammatory pathways, the short exposure (24 h), the microfluidic (physiological) context.

Figure D  Heatmap of top 6 DEGs across cell types and conditions
(NOT a biology figure, it is a transition figure between quantitative impact (Fig. 2A–C), and mechanistic, cell-type–specific biology (Figs. 3–6)).
So Fig. 2D should answer one question only: Are there shared transcriptional programs across immune cells, and how does cell-type specificity emerge on top of them?
Fig. 2D should contain ~25–30 genes total, split into three clearly interpretable blocks:
Block 1 — DEGs shared across multiple cell types (CORE)
This block visually justifies why Fig. 3 exists. Genes (12): RPL3, RPL7, RPL10, RPL13, RPLP0, RPS3, RPS6, RPS12, EIF3E, EIF4A1, NCBP1, NXF1. They are present across monocytes, B cells, CD4⁺, CD8⁺, NK; present across multiple conditions and diirectly linked to our strongest KEGG signal.
Block 2 — Antigen-handling genes (APC-biased but not exclusive). Genes (7): HLA-DRA, HLA-DRB1, CD74, CTSB, CTSD, LAMP1 and HLA-C. This prepares the reader for Fig. 4.
Block 3 — Minimal lineage markers (context, not mechanisms - not for interpretation, only to help read the heatmap). Genes: B cells: MS4A1, CD79A; T cells: CCR7, LTB, FCER2.
Structure:
Rows: genes (grouped by block, not clustered across blocks)
Columns: Cell type × condition
Values: Z-scored average expression (or log2FC)
Critical design rule: Cluster rows within blocks, not globally
You want readers to see:
shared ribosomal signal
APC-biased antigen genes
lineage context
Interpretation: “To bridge global DEG burden with gene-level responses, we visualized the expression of representative DEGs across immune cell types and exposure conditions (Fig. 2E). This analysis revealed a conserved transcriptional program dominated by translation and RNA-processing genes across immune populations, accompanied by more restricted modulation of antigen-handling genes in antigen-presenting cells and lineage-specific markers in adaptive immune cells.”

FIGURE 3. Conserved transcriptional programs induced by PSNP exposure

Figure 3 establishes the first biological principle of the paper: PSNP exposure induces a conserved transcriptional program across immune cell types, dominated by translation and RNA related processes.

This figure is:
program-centric
cross–cell type
non-mechanistic
non-toxicological


Fig 3A Dot plot: conserved translation / RNA program across immune cells
Purpose: to demonstrate existence and conservation of a shared transcriptional program across immune cell types and PSNP conditions.
This is the key panel of the entire manuscript.

Plot type: dot plot (Seurat-style or equivalent)
Y-axis: Genes
X-axis: Cell type × condition
Group by cell type
Within each group: Control | 40 nm | 200 nm | 40+200 nm
Encoding
Dot color: average log2 fold change (vs control)
Dot size: percentage of cells expressing the gene


Genes to include: translation / RNA processing:
RPL3
RPL7
RPL10
RPL13
RPS6
RPS12
EIF4A1
NCBP1
 (optional ninth: RPLP0 if visually needed)


These genes are present in DEG tables, appear across multiple cell types, support ribosome / RNA pathway enrichment, are mostly downregulated, consistent with Fig. 2C


What the reader should see
The same genes modulated across immune populations
Differences in magnitude, not presence/absence
No cell-type exclusivity

Fig 3B Heatmap: structure and modulation of the conserved program
Purpose: to show how the conserved translation/RNA program is modulated across cell types and PSNP sizes.This panel adds structure, not new biology.

Plot type: Heatmap
Axes
Rows: Genes (same genes as Fig. 3A, same order)
Columns: Cell type × condition


Values
Z-scored average expression
 (or Z-scored log2FC, but be consistent across figures)
Critical design rules
❗ Do NOT cluster genes globally
Keep gene order identical to Fig. 3A
You may cluster columns within cell types, but not required
Interpretation
Confirms coordinated modulation of translation-related genes
Highlights that the response is programmatic, not gene-by-gene noise


Fig 3C Representative volcano plot grounding the program
Purpose to demonstrate that the conserved program is visible in raw differential expression, not only after aggregation. (This is evidence, not discovery.)

Comparison
CD14⁺ monocytes: 40 nm vs control
(You may choose 200 nm instead if that is cleaner, but choose one.)
Plot type
Volcano plot
Axes
X-axis: log2 fold change
Y-axis: −log10(adjusted p-value)
Highlighting
Highlight only the same genes used in Fig. 3A/B:
RPL3, RPL7, RPL10, RPL13
RPS6, RPS12
EIF4A1, NCBP1
Use one color for highlighted genes.

What this panel proves: the conserved program is not an averaging artifact; these genes are among the strongest DEGs in at least one cell type



Figure 3 PSNP exposure induces a conserved translation- and RNA-related transcriptional program across immune cell types.
 (A) Dot plot showing the expression of representative translation and RNA-processing genes across immune cell types and PSNP exposure conditions. Dot size indicates the fraction of cells expressing each gene, and color represents average log2 fold change relative to control.
 (B) Heatmap of the same gene set, highlighting coordinated transcriptional modulation across cell types and exposure conditions.
 (C) Volcano plot of differential gene expression in CD14⁺ monocytes following exposure to 40 nm PSNPs, with translation- and RNA-related genes highlighted.

------------------------------------
SUPPLEMENTARY FIGURE S2  Analytical support for Figures 2 and 3
Supplementary Fig. S2 provides completeness and robustness for: DEG burden (Fig. 2), Directionality (Fig. 2C), Program identification (Fig. 3)


Fig. S2A Full volcano plots (supports: Fig. 2A, Fig. 2C, Fig. 3C)
What to show
Volcano plots for all cell types × all conditions
Same thresholds
Same axis limits
Minimal or no gene labels
Fig. S2B Extended DEG heatmaps (supports Fig. 2E, Fig. 3A/B) → samo za 40 ili 200 nm, tj. za isti condition koji je koriscen u Fig. 3. Moj predlog je da uradimo za oba, i onda izaberemo
What to show
Option A (preferred):
Heatmap of top 50 DEGs per major cell type
Monocytes
B cells
CD4⁺ T cells
Top 50 DEGs per cell type following 40 nm PSNP exposure
Structure:
Rows: genes (top 50 by adj. p or |logFC|)
Columns: cell type
No annotations
No pathway grouping
No interpretation.
Fig. S2C Full GO / KEGG enrichment results (supports Fig. 3 pathway interpretation)
What to show
Dot plots or bar plots of all enriched GO / KEGG terms
One panel per cell type
Include:
ribosome / RNA terms
metabolic terms
antigen handling terms
terms not discussed in main text


grouping KEGG terms by collapsed categories (same as main figures → exl table cleaned KEGG pathways), e.g.:
Translation / RNA
Metabolism / mitochondria
Antigen handling
Immune signaling (unemphasized)

Supplementary Fig. S2D — DEG overlap analysis (optional)
Supports Fig. 2 logic (shared vs condition-specific)
What to show: Venn diagrams (numbers only) One per major responding cell type (e.g. monocytes/B-cells). Omit entirely if overlap is trivial.


Supplementary Fig. S2 - legend
Supplementary Figure S2 | Differential expression and pathway analyses supporting Figures 2 and 3.
 (A) Volcano plots for all immune cell types and PSNP exposure conditions.
 (B) Extended heatmaps of differentially expressed genes per cell type.
 (C) Complete GO and KEGG enrichment results for each immune cell type.
 (D) Overlap analysis of differentially expressed genes across PSNP exposure conditions.

FIG 4  Size-dependent transcriptional remodeling in monocytes following PSNP exposure
Conceptual role: Fig. 4 provides the first cell-type–specific biological interpretation of the conserved transcriptional programs introduced in Figure 3.

Key message: Particle size modulates the magnitude and coherence of monocyte transcriptional responses, with 40 nm PSNPs inducing stronger and more coordinated transcriptional remodeling than 200 nm PSNPs.
Importantly: use words remodeling - not damage; coherence, not toxicity; stay within what DEGs + KEGG actually support
Fig. 4A Size-dependent differential gene expression in monocytes
What is shown
Two volcano plots (side by side):
Left: CD14⁺ monocytes, 40 nm vs control
Right: CD14⁺ monocytes, 200 nm vs control


Axes (identical for both plots)
X-axis: log2 fold change
Y-axis: −log10(adjusted p-value)
Genes highlighted (verified DEGs)
GSN (cytoskeletal organization)
MT-CO2 (mitochondrial metabolism)
TNFAIP3 (regulatory stress response)
SH3BGRL3 (regulatory / redox-associated; only if significant)


Interpretation (Results text language)
“Both particle sizes induced substantial transcriptional changes in monocytes; however, 40 nm PSNPs were associated with larger effect sizes and more pronounced modulation of genes involved in cytoskeletal organization, mitochondrial metabolism, and regulatory stress responses (Fig. 4A).”
Fig. 4B  Effect size distribution reveals stronger modulation by 40 nm PSNPs
Panel title
Figure 4B | Comparison of transcriptional effect sizes in monocytes
What is shown
Boxplot or violin plot of absolute log2 fold change values
Monocytes only
Axes
X-axis: PSNP condition (40 nm, 200 nm)
Y-axis: |log2FC|
Statistical annotation (optional but recommended)
Wilcoxon rank-sum test (p-value)
Interpretation
“Despite inducing a larger number of DEGs, 200 nm PSNPs elicited smaller transcriptional effect sizes compared to 40 nm particles, indicating that DEG burden does not directly reflect the magnitude of transcriptional remodeling (Fig. 4B).”
Fig. 4C Coherent pathway modulation in monocytes is size dependent
Panel title
Figure 4C | Size-dependent pathway coherence in monocytes
What is shown
Dot plot (preferred) or compact heatmap
Axes
X-axis: PSNP condition (40 nm, 200 nm, 40+200 nm)
Y-axis: Collapsed biological pathway themes
Collapsed pathway categories
Translation and RNA regulation
Cytoskeletal organization
Mitochondrial metabolic processes
Antigen processing and presentation
Cellular stress regulation
Encoding
Dot color: −log10(FDR)
Dot size: gene ratio
Interpretation
“Pathway enrichment analysis revealed that 40 nm PSNP exposure induced a more coherent modulation of multiple biological processes in monocytes, including cytoskeletal organization, mitochondrial metabolic processes, and antigen handling, whereas enrichment following 200 nm exposure was more diffuse and less coordinated (Fig. 4C).”
Fig. 4D Representative genes illustrate size-dependent transcriptional remodeling
Panel title
Figure 4D | Representative monocyte genes exhibit size-dependent transcriptional modulation
What is shown
Heatmap of selected genes across conditions
Axes
Rows: Genes
Columns: Control | 40 nm | 200 nm | 40+200 nm
Genes included
Cytoskeletal organization
GSN
Mitochondrial metabolism
MT-CO2
Antigen handling
CD74
HLA-DRA
Regulatory stress response
TNFAIP3
SH3BGRL3 (if significant)
Translation / RNA (contextual)
RPL7
RPS6
Design notes
Z-scored expression
Conditions not clustered
Gene clustering optional
Interpretation
“Representative genes associated with cytoskeletal organization, mitochondrial metabolism, antigen processing, and transcriptional regulation illustrate size-dependent transcriptional remodeling in monocytes across PSNP exposure conditions (Fig. 4D).”

Figure 4 | Size-dependent transcriptional remodeling in monocytes following PSNP exposure.
 (A) Volcano plots showing differential gene expression in CD14⁺ monocytes following exposure to 40 nm and 200 nm PSNPs. Selected genes associated with cytoskeletal organization, mitochondrial metabolism, and regulatory stress responses are highlighted.
 (B) Distribution of absolute log2 fold changes in monocytes, revealing larger transcriptional effect sizes following exposure to 40 nm PSNPs compared to 200 nm particles.
 (C) Pathway enrichment analysis highlighting size-dependent coherence in biological processes including translation and RNA regulation, cytoskeletal organization, mitochondrial metabolism, and antigen processing in monocytes.
 (D) Heatmap of representative monocyte genes illustrating size-dependent transcriptional modulation across PSNP exposure conditions.
This caption is final-submission ready.

SUPPLEMENTARY FIGURE S3 — Extended monocyte analyses supporting Figure 4
Conceptual role: Fig. S3 provides analytical completeness and transparency for the monocyte-specific analyses presented in Figure 4, without introducing additional biological interpretation.
Fig. S3A Full monocyte volcano plots
What is shown
Volcano plots for monocytes:
40 nm vs control
200 nm vs control
40+200 nm vs control
Purpose
Demonstrates full DEG distributions
Confirms robustness of Fig. 4A
Fig. S3B Complete pathway enrichment results in monocytes
What is shown
Full KEGG and/or GO enrichment results
One panel per condition
Purpose
Shows all enriched pathways, including those not discussed
Prevents cherry-picking concerns
Fig. S3C Gene-level effect size distributions at 40 and 200 nm PSNPs (ako se to vec vidi, onda je ovo nepotrebno)
What is shown
Boxplots (or violin plots) of log2FC for selected genes:
GSN (cytoskeletal organization)
MT-CO2 (mitochondrial metabolism)
TNFAIP3 (regulatory stress response)
CD74 (antigen processing)

Purpose: S3C shows gene-level effect-size distributions for a few representative monocyte genes, to support the claim that 40 nm PSNPs induce larger transcriptional modulation than 200 nm PSNPs.




Supplementary Figure S3 — FINAL CAPTION
Supplementary Figure S3 | Extended differential expression and pathway analyses in monocytes supporting Figure 4.
 (A) Volcano plots showing differential gene expression in monocytes across PSNP exposure conditions.
 (B) Complete KEGG and GO pathway enrichment results for monocytes.
 (C) Gene-level effect size distributions for selected monocyte genes.

FIG 5  Adaptive immune cells display restrained and lineage-specific transcriptional responses to PSNP exposure
Conceptual role of Figure 5 answers a new question, distinct from Figures 3 and 4: How do adaptive immune cells respond to PSNP exposure compared to innate immune cells?
 
Core message:
B cells and CD4⁺ T cells exhibit transcriptional responses that are more restrained and lineage-specific than those observed in monocytes, with limited pathway coherence and smaller effect sizes.
This is not a negative result — it is a biologically meaningful contrast.
Design principles for Figure 5
Fewer panels than Figure 4
Less mechanistic depth
Emphasis on specificity and restraint, not absence
Avoid inflammatory framing
Fig 5A DEG burden and effect size comparison in adaptive immune cells
Panel title
Figure 5A | Differential gene expression burden in adaptive immune cells
Grouped bar plot or dot plot
Axes
X-axis: Cell type (B cells, CD4⁺ T cells)
Y-axis: Number of DEGs (adjusted p < 0.05)
Encoding
Color-coded by condition:
40 nm
200 nm
40+200 nm
(Optional: outline bars or add markers for up/down split)
Purpose→ shows
Adaptive cells respond transcriptionally
Size dependence is present but muted
Interpretation (Results text)
“Notably, DEG burden alone did not reflect the degree of pathway-level organization, which differed markedly between adaptive and innate immune cells (see Fig. 6).”
Fig 5B Representative lineage-specific gene modulation (it is a qualitative, lineage-identity panel, not a “top DEG” panel.)
(for me: Its role is to show that:
adaptive immune cells retain their canonical identity
transcriptional changes are modest and lineage-specific
PSNP exposure does not induce a coherent reprogramming in these cells
This is why we do not choose: cytokines (IL1B, IL6, etc.), stress genes and mitochondrial genes. Those would pull the narrative toward toxicity or inflammation, which your pathway analysis does not support.)

Panel title
Figure 5B | Lineage-specific transcriptional modulation in adaptive immune cells
Plot type: Dot plot or compact heatmap
Axes
Y: Rows (genes)
X: Columns→ Cell type × condition
Encoding
Dot size = % expressing cells
Color = avg log2FC (relative to control)
Genes to include: 
B cells: 
MS4A1 (CD20): Canonical B-cell marker, Highly stable identity gene, Expressed broadly across B cells, Appears in DEG tables with modest modulation (not loss) → What it shows: B cells remain B cells after PSNP exposure. This is critical for arguing restraint.
CD79A: Core component of the B-cell receptor complex; Lineage-specific; Often slightly modulated but not suppressed. What it shows →Antigen-sensing machinery is preserved, not reprogrammed
CD74 (contextual, antigen handling). Links B cells to antigen presentation, Appears in KEGG “Antigen processing and presentation”, Already used in monocytes (Fig. 4), enabling comparison. Why optional → If CD74 modulation is weak in B cells, it can be dropped. Its role is contextual, not central
CD4⁺ T cells: 
CCR7: Why: Canonical naïve / central memory CD4⁺ T-cell marker; Cell-migration / homing gene; Sensitive to subtle perturbations. What it shows → T cells show modest modulation, not activation. Very useful to contrast with monocytes.
LTB: T-cell identity and immune architecture gene; Not an acute inflammatory cytokine; Stable expression in resting T cells. What it shows →CD4⁺ T cells maintain lineage programs
Shared transcriptional context
RPL7 and RPS6: Anchor Fig. 5 to Fig. 3 (translation/RNA program); Provide continuity across figures; Show that adaptive cells also participate in the shared program, but weakly. What they show → Shared transcriptional programs exist, but are not dominant in adaptive cells
What the reader should instantly see
Strong expression of lineage markers
Small shifts, not dramatic changes
No coherent stress/inflammatory signature
Purpose→ Shows:
Cell-type specificity
Preservation of lineage identity
Mild modulation rather than reprogramming
Interpretation
“Gene-level analysis revealed modest, lineage-specific transcriptional modulation in adaptive immune cells, with preservation of canonical B cell and CD4⁺ T cell identity markers (Fig. 5B).”

Fig 5C Pathway enrichment reveals limited coherence in adaptive immune cells
Panel title
Figure 5C | Pathway enrichment in adaptive immune cells following PSNP exposure

Plot type: Dot plot (preferred)
Axes
X-axis: Cell type (B cells, CD4⁺ T cells)
Y-axis: Collapsed pathway themes
Collapsed pathway themes. Use exactly the same categories as earlier, for consistency:
Translation and RNA regulation
Metabolic processes
Antigen processing and presentation
Cellular stress regulation
Encoding
Dot color: −log10(FDR)
Dot size: gene ratio
Purpose
Shows:
Enrichment exists
It is weaker and less coordinated than in monocytes
No dominant inflammatory program
Interpretation
“Pathway enrichment analysis in B cells and CD4⁺ T cells revealed modest modulation of transcriptional and metabolic processes, without the strong pathway coherence observed in monocytes (Fig. 5C).”

Figure 5 — FINAL CAPTION (LOCKED)
Figure 5 | Restrained and lineage-specific transcriptional responses in adaptive immune cells following PSNP exposure.
 (A) Differentially expressed gene counts in B cells and CD4⁺ T cells across PSNP exposure conditions.
 (B) Representative lineage-specific genes illustrating modest transcriptional modulation and preservation of adaptive immune cell identity.
 (C) Pathway enrichment analysis revealing limited and cell-type–specific transcriptional program modulation in adaptive immune cells.

SUPPLEMENTARY FIGURE S4 — Extended analyses of adaptive immune cells
Conceptual role: provides completeness and transparency for adaptive immune cell analyses supporting Figure 5.
Fig. S4A Full volcano plots
What is shown
Volcano plots for:
B cells (40, 200, mix)
CD4⁺ T cells (40, 200, mix)
Purpose
Shows full DEG distributions
Confirms restrained response
Fig. S4B Extended DEG heatmaps
Heatmaps of top 30–50 DEGs
y-axis: One per cell type (B cells, CD4⁺ T cells)
x-axis: One representative per every one of conditions (40, 200, mix)
Purpose
Demonstrates breadth without interpretation
Fig. S4C Complete pathway enrichment tables or plots
What is shown
All KEGG/GO enriched terms
One panel per cell type × condition
Purpose
Anti–cherry-picking
Supports Fig. 5C

Supplementary Figure S4 — FINAL CAPTION
Supplementary Figure S4 | Extended differential expression and pathway analyses in adaptive immune cells supporting Figure 5.
 (A) Volcano plots showing differential gene expression in B cells and CD4⁺ T cells across PSNP exposure conditions.
 (B) Extended heatmaps of differentially expressed genes in adaptive immune cells.
 (C) Complete pathway enrichment results for adaptive immune cell types.

FIGURE 6 — Integrative principles of immune transcriptional organization under complex PSNP exposure
Conceptual role of Figure 6: it answers the highest-level question of the paper: What organizing principles govern immune transcriptional responses to complex PSNP environments?
This figure must:
synthesize Figures 2–5
resolve apparent contradictions
elevate the work beyond descriptive scRNA-seq
Core take-home message (one sentence):
Immune responses to PSNP exposure are governed by transcriptional organization rather than DEG burden alone, and mixed particle exposure induces non-additive, emergent transcriptional states that cannot be inferred from single-size exposures.
This sentence is essentially our graphical abstract in words.
Fig 6A DEG burden does not predict pathway organization
Panel title
Figure 6A | Differential gene expression burden does not linearly predict pathway enrichment

Plot type:Scatter plot
Axes
X-axis: Number of DEGs
Y-axis: Number of significantly enriched KEGG pathways
Data points
Each dot = one cell type × one condition
Include: (ajde stavi sve za pocetak, pa cemo izabrati samo dva celijska tipa, ako je panel prekomplikovan, mozda monocyte i B ili CD4+...)
Monocytes
B cells
CD4⁺ T cells
(Optionally CD8⁺ T cells)
Encoding
Color: cell type
Shape: exposure condition
○ 40 nm
△ 200 nm
◇ 40 + 200 nm
What the reader should immediately see
Monocytes:
 → modest DEG numbers
 → high pathway counts
Adaptive cells:
 → high DEG numbers (in some conditions)
 → low pathway counts
Mixed exposure deviates from linear trends
Interpretation:
“Across immune cell types, the number of differentially expressed genes did not correlate linearly with the number of enriched biological pathways, indicating fundamental differences in transcriptional organization between innate and adaptive immune responses (Fig. 6A).”
Fig 6B Mixed exposure induces non-additive transcriptional responses
Panel title
Figure 6B | Combined PSNP exposure elicits non-additive transcriptional responses
Plot type: Paired dot plot or grouped bar plot
What is compared
For each selected cell type (recommend Monocytes + B cells):
Expected additive response
 = mean of (40 nm + 200 nm)
Observed mixed response
 = 40 + 200 nm
Metric: Number of enriched pathways (strongest, most conceptual)
Axes
X-axis: Expected vs observed
Y-axis: Pathway count
What the reader should see
Mixed exposure ≠ expected additive response
Often attenuated or qualitatively different
Clear deviation from linearity
Interpretation (Results-ready wording)
“The transcriptional response to combined 40 + 200 nm PSNP exposure deviated from the expected additive effect of individual particle sizes, revealing non-linear and emergent transcriptional states across immune cell types (Fig. 6B).”
Fig 6C  Conceptual model of immune transcriptional organization under PSNP exposure
Panel title
Figure 6C | Conceptual model of immune transcriptional organization under single-size and mixed PSNP exposure
Plot type: Schematic / illustration (no data)

What it should depict
Three exposure scenarios:
40 nm
Strong, coherent transcriptional remodeling
Innate-biased
Program-level organization
200 nm
Broad but diffuse transcriptional perturbation
Higher DEG dispersion
Lower pathway coherence
40 + 200 nm
Emergent, non-additive state
Attenuated or selective program activation
Not predictable from single-size exposures


FIGURE 6 — FINAL CAPTION (LOCKABLE)
Figure 6 | Integrative principles of immune transcriptional organization under complex PSNP exposure.
 (A) Relationship between differential gene expression burden and pathway enrichment across immune cell types and PSNP exposure conditions, demonstrating that DEG number does not linearly predict biological pathway organization.
 (B) Comparison of expected additive and observed transcriptional responses following combined 40 + 200 nm PSNP exposure, revealing non-additive and emergent immune transcriptional states.
 (C) Conceptual model summarizing size-dependent and non-additive immune transcriptional responses to single-size and mixed PSNP exposure.

SUPPLEMENTARY FIGURE S5 — Supporting analyses for Figure 6
Conceptual role: provides quantitative depth and robustness for Figure 6 without cluttering the main narrative.
Fig. S5A Expanded DEG vs pathway scatter
What it shows
Same as Fig. 6A
Includes:
all cell types
all conditions
regression lines (optional)
Fig. S5B Non-additivity metrics
What it shows
Expected vs observed comparisons using: number of enriched KEGG pathways or DEG counts 
How to plot it: Paired dot plot or grouped bar plot
For each selected cell type:
Condition
Value
Expected additive
(40 nm + 200 nm) / 2
Observed mix
40 + 200 nm

So visually:
Expected (open symbol or light bar)
Observed (filled symbol or dark bar)

Supports Fig. 6B robustness.

Fig. S5C Cell-type–resolved non-additivity (Treba da prikazemo sve celijske tipove i sve conditions)
What it shows
Separate panels per cell type
Shows that non-additivity is not monocyte-exclusive



Supplementary Figure S5 — Caption
Supplementary Figure S5 | Extended analyses supporting integrative transcriptional principles shown in Figure 6.

Razlika izmedju 6B, S5B i S5C
Start from the anchor: Figure 6B
What Fig. 6B does
Uses ONE metric (e.g. number of enriched KEGG pathways)
Uses 1–2 representative cell types (e.g. monocytes + B cells)
Shows:
expected additive response
observed mixed response
Makes the core claim: Mixed exposure is non-additive
Fig. 6B is the clean, illustrative example.

Everything in Supplementary Fig. S5 exists only to defend Fig. 6B.
Fig. S5B 
What question S5B answers “Is non-additivity dependent on the metric we chose?”
What if someone says: “Your non-additivity is just because you chose pathway count”?
What S5B does
Keeps the same cell types as Fig. 6B (monocytes + B cells)
Keeps expected vs observed logic
Changes ONLY the metric
Example:
Fig. 6B → pathway count
S5B → mean |log2FC| (or DEG count, if you prefer)
Interpretation → If non-additivity is still visible: It is not an artifact of the chosen readout. So S5B tests metric robustness
Fig. S5C What question S5C answers“Is non-additivity specific to one cell type?”
What if someone says: “This is true for monocytes, but not for other immune cells”?
What S5C does
Keeps the SAME metric as Fig. 6B
 (e.g. pathway count)
Keeps expected vs observed logic
Changes ONLY the cell type
So you show:
Monocytes
B cells
CD4⁺ T cells
(optionally CD8⁺)
Each cell type gets its own mini comparison:
expected additive
observed mixed
Interpretation → If deviations are present in more than one cell type: Non-additivity is system-level, not cell-specific. So S5C tests biological generality

Detailed explanation:
Fig. S5C is a small multi-panel figure, for example:
S5C
 ├── Monocytes
 ├── B cells
 ├── CD4+ T cells
 └── (optional) CD8+ T cells

Each mini-panel is identical in structure, only the cell type changes.
What each mini-panel contains
For each cell type, show two values only:
Expected additive response
Calculated as:
 (40 nm + 200 nm) / 2
Using the same metric as Fig. 6B
 (e.g. number of enriched KEGG pathways)
Observed mixed response
Value from 40 + 200 nm condition
How to plot it
Option A (simplest, recommended): Paired bar plot
For each cell type:
Light bar = expected additive
Dark bar = observed mixed
Y-axis: Pathway count (or whatever metric you used in Fig. 6B)
X-axis: Expected vs Observed
Option B (even cleaner): Paired dot plot
For each cell type:
Open circle = expected
Filled circle = observed
Connected by a line
This visually emphasizes deviation from expectation.

Side-by-side comparison (this usually makes it click)
Figure
What changes?
What stays fixed?
What it proves
Fig. 6B
nothing
one metric, 1–2 cell types
main claim
S5B
metric
same cell types
robustness to measurement
S5C
cell type
same metric
robustness across biology
