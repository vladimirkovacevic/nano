# make_enrichment_summary.py
# Pravi sumarni Excel/CSV i 3 figure na osnovu:
#  - g:Profiler rezultata (XLSX)
#  - CSV sa genima i logFC (STRING/Cytoscape izlaz)

import os, math, re
import pandas as pd
import matplotlib.pyplot as plt

# ========== PUTANJE (PRILAGOĐENO TVOM RASPOREDU) ==========
ENR_SUMMARY_DIR = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Pathway Enrichment\Enrichment Summary"

# g:Profiler rezultat – TAČNA LOKACIJA koju si naveo:
ENR_FILE = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Pathway Enrichment\Bcell_genes_gprofiler_results.xlsx"

# Ulazni CSV sa genima/logFC (iz Cytoscape/STRING koraka):
INPUT_CSV = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Cytoscape\B_cell_sample_1_ImmuneMeta_overlap_STRING.csv"

# Izlazi:
OUT_XLSX   = os.path.join(ENR_SUMMARY_DIR, "Bcell_genes_enrichment_summary.xlsx")
OUT_CSV    = os.path.join(ENR_SUMMARY_DIR, "Bcell_genes_enrichment_summary.csv")
FIG_BAR     = os.path.join(ENR_SUMMARY_DIR, "FigA_Bcell_bar_logFC.png")
FIG_BUBBLE  = os.path.join(ENR_SUMMARY_DIR, "FigB_Bcell_bubble_Top5.png")
FIG_NETWORK = os.path.join(ENR_SUMMARY_DIR, "FigC_Bcell_term_gene_network_Top5.png")

# Nazivi kolona koje očekujemo u INPUT_CSV:
GENE_COL  = "gene_name"
LOGFC_COL = "logfoldchanges"
TOP_N     = 5

# ========== POMOĆNE ==========
def sniff_sep(path):
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        head = "".join([next(f, "") for _ in range(3)])
    cands = [",", "\t", ";", "|"]
    counts = {sep: head.count(sep) for sep in cands}
    sep = max(counts, key=counts.get)
    return sep if counts[sep] > 0 else ","

def neglog10(p):
    try:
        if p is None or p <= 0: return 0.0
        return -math.log10(p)
    except Exception:
        return 0.0

def classify_term(name, src):
    n = str(name).lower()
    s = str(src).upper()
    if any(k in n for k in ["interferon", "antiviral", "viral", "ifn"]):
        return "Interferon/antiviral"
    if any(k in n for k in ["tnf", "nf-kb", "nfkb", "toll-like", "tlr", "jak-stat", "cytokine", "chemokine", "il-"]):
        return "Immune signaling"
    if any(k in n for k in ["antigen", "mhc", "hla", "presentation"]):
        return "Antigen presentation"
    if any(k in n for k in ["apoptosis", "cell death", "pyroptosis", "necroptosis", "autophagy"]):
        return "Cell death"
    if s == "REAC" and any(k in n for k in ["metabolism", "metabolic", "biosynthesis", "catabolism"]):
        return "Metabolic"
    return "Other"

def parse_intersections(x):
    if isinstance(x, (list, tuple, set)):
        return [str(g).strip() for g in x if str(g).strip()]
    if isinstance(x, str):
        parts = x.replace(";", ",").split(",")
        return [t.strip() for t in parts if t.strip()]
    return []

def pick(cols, cands):
    for c in cands:
        if c in cols: return c
    # tolerantno na razmake/veličinu slova
    norm = {re.sub(r"\s+", "", c).lower(): c for c in cols}
    for c in cands:
        key = re.sub(r"\s+", "", c).lower()
        if key in norm: return norm[key]
    return None

# ========== PROVERE ==========
os.makedirs(ENR_SUMMARY_DIR, exist_ok=True)
if not os.path.exists(ENR_FILE):
    raise SystemExit(f"[ERR] Ne postoji ENR_FILE:\n{ENR_FILE}")
if not os.path.exists(INPUT_CSV):
    raise SystemExit(f"[ERR] Ne postoji INPUT_CSV:\n{INPUT_CSV}")

# ========== 1) UČITAJ GENE + logFC ==========
sep = sniff_sep(INPUT_CSV)
df_in = pd.read_csv(INPUT_CSV, sep=sep, encoding="utf-8-sig")
if GENE_COL not in df_in.columns:
    raise SystemExit(f"[ERR] Kolona '{GENE_COL}' ne postoji u {INPUT_CSV}. Kolone: {list(df_in.columns)}")
if LOGFC_COL not in df_in.columns:
    # fallback ako je drugačije ime kolone
    alt = None
    for c in ["log2fc", "logfc"]:
        if c in df_in.columns: alt = c; break
    if alt is None:
        raise SystemExit(f"[ERR] Kolona '{LOGFC_COL}' ne postoji u {INPUT_CSV}. Kolone: {list(df_in.columns)}")
    LOGFC_COL = alt

df_in["_gene_upper"] = df_in[GENE_COL].astype(str).str.strip().str.upper()
logfc_map = dict(zip(df_in["_gene_upper"], df_in[LOGFC_COL]))

# ========== 2) UČITAJ g:Profiler (sheet 'Enrichment_all' ili prvi) ==========
try:
    res = pd.read_excel(ENR_FILE, sheet_name="Enrichment_all")
except ValueError:
    res = pd.read_excel(ENR_FILE, sheet_name=0)

cols = list(res.columns)
term_id_col   = pick(cols, ["Term_ID","term_id","native","term_id_in_database","id"])
term_name_col = pick(cols, ["Term","term_name","name","term"])
src_col       = pick(cols, ["Source","source"])
padj_col      = pick(cols, ["adj_p_value","adjusted_p_value","padj","qvalue","q_value","adjusted_p_value_exclusive_fdr"])
intrs_col     = pick(cols, ["intersections","intersection","intersect","genes","intersections_found"])
hits_col      = pick(cols, ["hits","intersection_size"])
neglog_col    = pick(cols, ["-log10(adj_p)","negative_log10_of_adjusted_p_value","neglog_adj"])
p_col         = pick(cols, ["p_value","p-value","pvalue"])

rename = {}
if term_id_col:   rename[term_id_col]   = "Term_ID"
if term_name_col: rename[term_name_col] = "Term"
if src_col:       rename[src_col]       = "Source"
if padj_col:      rename[padj_col]      = "adj_p_value"
if intrs_col:     rename[intrs_col]     = "intersections"
if hits_col:      rename[hits_col]      = "hits"
if neglog_col:    rename[neglog_col]    = "-log10(adj_p)"
if p_col:         rename[p_col]         = "p_value"
if rename:
    res = res.rename(columns=rename)

# Fallback kolone
if "adj_p_value" not in res.columns:
    if "-log10(adj_p)" in res.columns:
        res["adj_p_value"] = 10 ** (-res["-log10(adj_p)"])
    else:
        res["adj_p_value"] = res.get("p_value", pd.Series([None]*len(res)))
if "Term_ID" not in res.columns:
    res["Term_ID"] = (res.get("Source","SRC").astype(str) + ":" + res.get("Term","TERM").astype(str))
if "hits" not in res.columns:
    res["hits"] = res.get("intersections","").apply(parse_intersections).apply(len)
if "-log10(adj_p)" not in res.columns:
    res["-log10(adj_p)"] = res["adj_p_value"].apply(neglog10)
if "Source" not in res.columns: res["Source"] = "NA"
if "Term" not in res.columns:   res["Term"]   = res["Term_ID"]

# Kategorije i rangiranje
res["Category"] = [classify_term(t, s) for t, s in zip(res["Term"], res["Source"])]
res_sorted = res.sort_values(["adj_p_value","p_value"], ascending=[True, True], na_position="last").reset_index(drop=True)

# ========== 3) Top N + membership ==========
top_terms = res_sorted.head(TOP_N).copy()

summary_rows, membership_rows = [], []
for idx, r in top_terms.iterrows():
    term_id = r["Term_ID"]; term = r["Term"]; source = r["Source"]; cat = r["Category"]
    adjp = r["adj_p_value"]; hits = int(r.get("hits", 0))
    genes_in_term = parse_intersections(r.get("intersections", ""))

    genes_dir = []
    for g in genes_in_term:
        g_up = str(g).upper().strip()
        fc = logfc_map.get(g_up, float("nan"))
        if pd.isna(fc):
            direction = "NA"; label = f"{g}"
        else:
            direction = "↑" if fc > 0 else "↓" if fc < 0 else "0"
            label = f"{g} ({direction} {fc:.2f})"
        genes_dir.append(label)
        membership_rows.append({
            "Term_ID": term_id, "Term": term, "Source": source, "Category": cat,
            "Gene": g, "logFC": None if pd.isna(fc) else fc, "Direction": direction
        })

    summary_rows.append({
        "Rank": idx + 1, "Term_ID": term_id, "Term": term, "Source": source, "Category": cat,
        "adj_p_value": adjp, "-log10(adj_p)": neglog10(adjp), "hits": hits,
        "Genes (with direction & logFC)": "; ".join(genes_dir) if genes_dir else ""
    })

summary_df  = pd.DataFrame(summary_rows)
member_long = pd.DataFrame(membership_rows)

# ========== 4) Snimi izlaze ==========
with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as wr:
    summary_df.to_excel(wr, index=False, sheet_name="Summary_Top5")
    member_long.to_excel(wr, index=False, sheet_name="Membership_Long")
summary_df.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

print("[OK] Sačuvano:")
print("  ", OUT_XLSX)
print("  ", OUT_CSV)

# ========== 5) Figure ==========
# A) Bar plot logFC (svi jedinstveni geni)
genes_fc = (df_in[[GENE_COL, LOGFC_COL]]
            .dropna()
            .drop_duplicates(subset=[GENE_COL])
            .sort_values(LOGFC_COL))

plt.figure(figsize=(8, 6))
plt.barh(genes_fc[GENE_COL], genes_fc[LOGFC_COL])
plt.axvline(0, linestyle="--")
plt.xlabel("log2 fold change"); plt.ylabel("Gene")
plt.title("B cells: log2 fold changes")
plt.tight_layout(); plt.savefig(FIG_BAR, dpi=300); plt.close()
print("[OK] Snimljeno:", FIG_BAR)

# B) Bubble plot Top5 termina
plot_df = top_terms.copy()
plot_df["-log10(adj_p)"] = plot_df["adj_p_value"].apply(neglog10)
plt.figure(figsize=(10, 6))
y = plot_df["Term"][::-1]
x = plot_df["-log10(adj_p)"][::-1]
s = (plot_df["hits"][::-1] * 160).clip(lower=80)
plt.scatter(x, y, s=s)
plt.xlabel("-log10(adj p)"); plt.ylabel("Term")
plt.title("Top 5 enriched terms (FDR)")
plt.tight_layout(); plt.savefig(FIG_BUBBLE, dpi=300); plt.close()
print("[OK] Snimljeno:", FIG_BUBBLE)

# C) Term–Gene bipartite (Top5)
edges, genes_set = [], set()
for _, r in top_terms.iterrows():
    term = str(r["Term"])
    for g in parse_intersections(r.get("intersections", "")):
        edges.append((term, g)); genes_set.add(g)

terms = list(dict.fromkeys([e[0] for e in edges]))
genes_unique = sorted(list(genes_set))

def linspace(n):
    if n <= 1: return [0.5]
    step = 1.0 / (n - 1); return [1.0 - i*step for i in range(n)]

term_y = linspace(len(terms)); gene_y = linspace(len(genes_unique))
term_pos = {t: (0.0, term_y[i]) for i, t in enumerate(terms)}
gene_pos = {g: (1.0, gene_y[i]) for i, g in enumerate(genes_unique)}

plt.figure(figsize=(11, 6))
for t, g in edges:
    x1, y1 = term_pos[t]; x2, y2 = gene_pos[g]
    plt.plot([x1, x2], [y1, y2])
for t, (x, y) in term_pos.items():
    plt.scatter([x], [y], marker="s", s=220); plt.text(x-0.02, y, t, ha="right", va="center")
for g, (x, y) in gene_pos.items():
    plt.scatter([x], [y], s=90); plt.text(x+0.02, y, g, ha="left", va="center")

plt.axis("off"); plt.title("Term–Gene bipartite network (Top 5 terms)")
plt.tight_layout(); plt.savefig(FIG_NETWORK, dpi=300); plt.close()
print("[OK] Snimljeno:", FIG_NETWORK)
