# gprofiler_enrichment_allinone.py  (robustne kolone)

import os, math
import pandas as pd
import matplotlib.pyplot as plt
from gprofiler import GProfiler

IN_CSV  = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Pathway Enrichment\B_cell_sample_1_ImmuneMeta_overlap_STRING.csv"
OUT_DIR = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Pathway Enrichment"

GENE_COL  = "gene_name"
LOGFC_COL = "logfoldchanges"
SOURCES   = ["REAC", "WP", "GO:BP", "GO:MF", "GO:CC", "DO", "HP"]

def neglog10(p):
    try:
        if p is None or p <= 0: return 0.0
        import math
        return -math.log10(p)
    except Exception:
        return 0.0

def classify_term(row):
    name = str(row.get("Term", "")).lower()
    src  = str(row.get("Source", "")).upper()
    if any(k in name for k in ["interferon", "antiviral", "viral", "ifn"]):
        return "Interferon/antiviral"
    if any(k in name for k in ["tnf", "nf-kb", "nfkb", "toll-like", "tlr", "jak-stat", "cytokine", "chemokine", "il-"]):
        return "Immune signaling"
    if any(k in name for k in ["antigen", "mhc", "hla", "presentation"]):
        return "Antigen presentation"
    if any(k in name for k in ["apoptosis", "cell death", "pyroptosis", "necroptosis", "autophagy"]):
        return "Cell death"
    if src == "REAC" and any(k in name for k in ["metabolism", "metabolic", "biosynthesis", "catabolism"]):
        return "Metabolic"
    return "Other"

def count_hits(x):
    if isinstance(x, (list, tuple, set)): return len(x)
    if isinstance(x, str): return len([t for t in x.split(",") if t.strip()])
    return 1

# --- 1) Učitaj gene ---
os.makedirs(OUT_DIR, exist_ok=True)
if not os.path.exists(IN_CSV):
    raise SystemExit(f"[ERR] Ne postoji ulazni CSV: {IN_CSV}")

df_in = pd.read_csv(IN_CSV)
if GENE_COL not in df_in.columns:
    raise SystemExit(f"[ERR] Kolona '{GENE_COL}' ne postoji. Kolone: {list(df_in.columns)}")

genes = (df_in[GENE_COL].astype(str).str.strip()
         .replace({"": pd.NA}).dropna().unique().tolist())
if not genes:
    raise SystemExit("[ERR] Ulazna lista gena je prazna.")
print(f"[INFO] Ulaznih gena: {len(genes)} → {genes}")

# --- 2) g:Profiler ---
gp = GProfiler(return_dataframe=True)
res = gp.profile(
    organism="hsapiens",
    query=genes,
    sources=SOURCES,
    user_threshold=1.0,
    all_results=True,
    no_evidences=False
)
if res is None or res.empty:
    raise SystemExit("[WARN] g:Profiler vratio prazan skup.")

# --- 2a) Normalizuj razne varijante kolona ---
# Kandidati za ID, ime termina, p-vrednosti, intersections
TERM_ID_CAND = ["Term_ID","term_id","native","term_id_in_database","id"]
TERM_NAME_CAND = ["Term","term_name","name","term"]
SRC_CAND = ["Source","source"]
PV_CAND = ["p_value","p-value","pvalue"]
ADJPV_CAND = ["adj_p_value","adjusted_p_value","padj","qvalue","q_value"]
NEGLOG_ADJ_CAND = ["negative_log10_of_adjusted_p_value","-log10(adj_p)"]
INTRS_CAND = ["intersections","intersection","intersect","genes"]

def pick(cols, cands):
    for c in cands:
        if c in cols: return c
    return None

cols = list(res.columns)

term_id_col   = pick(cols, TERM_ID_CAND)
term_name_col = pick(cols, TERM_NAME_CAND)
src_col       = pick(cols, SRC_CAND)
p_col         = pick(cols, PV_CAND)
padj_col      = pick(cols, ADJPV_CAND)
neglog_col    = pick(cols, NEGLOG_ADJ_CAND)
intrs_col     = pick(cols, INTRS_CAND)

# Preimenuj u standard
rename_map = {}
if term_id_col:   rename_map[term_id_col]   = "Term_ID"
if term_name_col: rename_map[term_name_col] = "Term"
if src_col:       rename_map[src_col]       = "Source"
if p_col:         rename_map[p_col]         = "p_value"
if padj_col:      rename_map[padj_col]      = "adj_p_value"
if neglog_col:    rename_map[neglog_col]    = "neglog_adj"
if intrs_col:     rename_map[intrs_col]     = "intersections"
if rename_map:
    res = res.rename(columns=rename_map)

# Ako nema adj_p_value, rekonstruiši iz neglog_adj ako postoji
if "adj_p_value" not in res.columns:
    if "neglog_adj" in res.columns:
        res["adj_p_value"] = 10 ** (-res["neglog_adj"])
    else:
        res["adj_p_value"] = res.get("p_value", pd.Series([None]*len(res)))

# Ako nema Term_ID, napravi ga iz "Term" + "Source" kao fall-back
if "Term_ID" not in res.columns:
    res["Term_ID"] = (res.get("Source", "SRC") + ":" + res.get("Term", "TERM").astype(str))

# Derivati
res["-log10(adj_p)"] = res["adj_p_value"].apply(neglog10)
if "intersections" not in res.columns:
    res["intersections"] = ""
res["hits"] = res["intersections"].apply(count_hits)
if "Source" not in res.columns:
    res["Source"] = "NA"
if "Term" not in res.columns:
    res["Term"] = res.get("Term_ID")

res["Category"] = res.apply(classify_term, axis=1)

# --- 3) Sačuvaj tabele ---
res_sorted = res.sort_values(["adj_p_value", "p_value"], ascending=[True, True], na_position="last").reset_index(drop=True)

csv_all  = os.path.join(OUT_DIR, "Bcell_genes_gprofiler_all.csv")
xlsx_out = os.path.join(OUT_DIR, "Bcell_genes_gprofiler_results.xlsx")

res_sorted.to_csv(csv_all, index=False, encoding="utf-8-sig")
with pd.ExcelWriter(xlsx_out, engine="openpyxl") as wr:
    cols_in = [GENE_COL] + ([LOGFC_COL] if LOGFC_COL in df_in.columns else [])
    df_in[cols_in].to_excel(wr, index=False, sheet_name="Input_genes")
    res_sorted.to_excel(wr, index=False, sheet_name="Enrichment_all")
    # Top 10 po izvoru
    tops = []
    for src in SOURCES:
        sub = res_sorted[res_sorted["Source"]==src].head(10).copy()
        sub.insert(0, "Rank_in_source", range(1, len(sub)+1))
        tops.append(sub)
    if tops:
        pd.concat(tops, ignore_index=True).to_excel(wr, index=False, sheet_name="Top_by_source")
    # Gene↔Term matrica (Top 30 termina)
    top_terms = res_sorted.head(30)["Term_ID"].tolist()
    terms_df = res_sorted[res_sorted["Term_ID"].isin(top_terms)][["Term_ID","Term","Source","intersections"]].copy()
    rows = []
    for _, r in terms_df.iterrows():
        ints = r["intersections"]
        if isinstance(ints, str):
            glist = [g.strip() for g in ints.split(",") if g.strip()]
        elif isinstance(ints, (list, tuple, set)):
            glist = list(ints)
        else:
            glist = []
        for g in glist:
            rows.append({"Term_ID":r["Term_ID"], "Term":r["Term"], "Source":r["Source"], "Gene":g})
    mem_df = pd.DataFrame(rows)
    if not mem_df.empty:
        mem_df["val"] = 1
        mat = mem_df.pivot_table(index=["Term_ID","Term","Source"], columns="Gene", values="val", fill_value=0)
        mat.to_excel(wr, sheet_name="GeneTerm_matrix")

print(f"[OK] Sačuvano:\n  {csv_all}\n  {xlsx_out}")

# --- 4) Bubble plot (Top 15 po FDR) ---
dbg_res = os.path.join(OUT_DIR, "DEBUG_gprofiler_raw_columns.csv")
res_sorted.head(50).to_csv(dbg_res, index=False, encoding="utf-8-sig")

plot_df = res_sorted[res_sorted["adj_p_value"].notna()].head(15).copy()
if not plot_df.empty:
    plt.figure(figsize=(10, 6))
    y = plot_df["Term"][::-1]
    x = plot_df["-log10(adj_p)"][::-1]
    s = (plot_df["hits"][::-1] * 120).clip(lower=60)
    plt.scatter(x, y, s=s)
    plt.xlabel("-log10(adj p)")
    plt.ylabel("Term")
    plt.title("Pathway / Process enrichment – Top terms (B cell, genes)")
    plt.tight_layout()
    png_out = os.path.join(OUT_DIR, "Bcell_genes_gprofiler_bubble.png")
    plt.savefig(png_out, dpi=300)
    plt.close()
    print(f"[OK] Bubble plot: {png_out}")
else:
    print("[WARN] Nema dovoljno termina za bubble plot.")

# --- 5) Cytoscape mreže ---
# (a) Term↔Gene edges (za SVE termine)
rows = []
for _, r in res_sorted.iterrows():
    ints = r.get("intersections", [])
    if isinstance(ints, str):
        glist = [g.strip() for g in ints.split(",") if g.strip()]
    elif isinstance(ints, (list, tuple, set)):
        glist = list(ints)
    else:
        glist = []
    for g in glist:
        rows.append({"Term_ID": r["Term_ID"], "Term": r["Term"], "Source": r["Source"], "Gene": g})

mem_all = pd.DataFrame(rows)
edges_out = os.path.join(OUT_DIR, "Bcell_genes_TermGene_edges.tsv")
if not mem_all.empty:
    mem_all.to_csv(edges_out, sep="\t", index=False)
    print(f"[OK] Cytoscape edges (Term-Gene): {edges_out}")
else:
    print("[WARN] Nije napravljen Term–Gene edges TSV.")

# (b) Term↔Term overlap (Jaccard) za termine sa ≥2 gena
def jaccard(a,b):
    a, b = set(a), set(b)
    if not a and not b: return 0.0
    return len(a & b) / max(1, len(a | b))

term2genes = {}
for _, r in mem_all.iterrows():
    term2genes.setdefault(r["Term_ID"], set()).add(r["Gene"])

valid_tids = [tid for tid, gl in term2genes.items() if len(gl) >= 2]
rows_tt = []
vt = list(valid_tids)
for i in range(len(vt)):
    for j in range(i+1, len(vt)):
        t1, t2 = vt[i], vt[j]
        w = jaccard(term2genes[t1], term2genes[t2])
        if w > 0:
            rows_tt.append({"Term_ID_A": t1, "Term_ID_B": t2, "Jaccard": w})

tt_df = pd.DataFrame(rows_tt)
if not tt_df.empty:
    id2name = dict(zip(res_sorted["Term_ID"], res_sorted["Term"]))
    tt_df["Term_A"] = tt_df["Term_ID_A"].map(id2name)
    tt_df["Term_B"] = tt_df["Term_ID_B"].map(id2name)
    tt_out = os.path.join(OUT_DIR, "Bcell_genes_TermTerm_overlap.tsv")
    tt_df.to_csv(tt_out, sep="\t", index=False)
    print(f"[OK] Cytoscape Term–Term overlap: {tt_out}")
else:
    print("[WARN] Nije napravljen Term–Term overlap TSV (premalo preklapanja ili 1-gen termini).")
