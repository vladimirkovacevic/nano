# make_publication_text_summary.py
import os, re, math
import pandas as pd

# === PODESI PUTANJE (po tvom rasporedu) ===
ENR_XLSX  = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Pathway Enrichment\Bcell_genes_gprofiler_results.xlsx"
INPUT_CSV = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Cytoscape\B_cell_sample_1_ImmuneMeta_overlap_STRING.csv"
OUT_TXT   = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Pathway Enrichment\Enrichment Summary\Bcell_genes_publication_summary.txt"

DATASET_NAME = "B cell sample 1"
FDR_CUTOFF   = 0.05
TOP_N        = 5

def sniff_sep(path):
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        head = "".join([next(f, "") for _ in range(3)])
    cands = [",", "\t", ";", "|"]
    counts = {sep: head.count(sep) for sep in cands}
    return max(counts, key=counts.get) if max(counts.values())>0 else ","

def pick(cols, cands):
    for c in cands:
        if c in cols: return c
    norm = {re.sub(r"\s+", "", c).lower(): c for c in cols}
    for c in cands:
        key = re.sub(r"\s+", "", c).lower()
        if key in norm: return norm[key]
    return None

def parse_intersections(x):
    if isinstance(x, str):
        parts = x.replace(";", ",").split(",")
        return [t.strip() for t in parts if t.strip()]
    return []

def neglog10(p):
    try:
        return -math.log10(p) if p and p>0 else 0.0
    except:
        return 0.0

# --- učitaj gene + logFC ---
sep = sniff_sep(INPUT_CSV)
genes = pd.read_csv(INPUT_CSV, sep=sep, encoding="utf-8-sig")
gene_col  = pick(genes.columns, ["gene_name","gene","symbol"])
logfc_col = pick(genes.columns, ["logfoldchanges","log2fc","logfc"])
if gene_col is None or logfc_col is None:
    raise SystemExit("[ERR] Nema gene/logFC kolona u INPUT_CSV.")

genes["_UP"] = genes[gene_col].astype(str).str.strip().str.upper()
logfc_map = dict(zip(genes["_UP"], genes[logfc_col]))

# --- učitaj g:Profiler rezultate ---
try:
    enr = pd.read_excel(ENR_XLSX, sheet_name="Enrichment_all")
except ValueError:
    enr = pd.read_excel(ENR_XLSX, sheet_name=0)

# mapiraj kolone
cols = list(enr.columns)
term_id   = pick(cols, ["Term_ID","term_id","native","id"])
term_name = pick(cols, ["Term","term_name","name","term"])
source    = pick(cols, ["Source","source"])
padj      = pick(cols, ["adj_p_value","adjusted_p_value","padj","qvalue","q_value"])
pval      = pick(cols, ["p_value","pvalue","p-value"])
hits_col  = pick(cols, ["hits","intersection_size"])
intrs     = pick(cols, ["intersections","intersection","genes","intersections_found"])

rename = {}
if term_id:   rename[term_id]   = "Term_ID"
if term_name: rename[term_name] = "Term"
if source:    rename[source]    = "Source"
if padj:      rename[padj]      = "adj_p_value"
if pval:      rename[pval]      = "p_value"
if hits_col:  rename[hits_col]  = "hits"
if intrs:     rename[intrs]     = "intersections"
if rename:    enr = enr.rename(columns=rename)

for need in ["Term","Source","adj_p_value","intersections"]:
    if need not in enr.columns:
        raise SystemExit(f"[ERR] U Enrichment_all fali kolona: {need}")

enr["hits"] = enr.get("hits", enr["intersections"].apply(parse_intersections).apply(len))
enr = enr.sort_values(["adj_p_value","p_value"], ascending=[True, True], na_position="last").reset_index(drop=True)

# --- filtriraj & izdvoj top N ---
sig = enr[enr["adj_p_value"] <= FDR_CUTOFF].copy()
n_sig = len(sig)
top = (sig if n_sig>0 else enr).head(TOP_N).copy()

# --- sastavi neutralni tekst ---
lines = []
lines.append(f"{DATASET_NAME} — funkcionalno obogaćenje (automatizovani rezime)")
lines.append("")
lines.append("Ulazni podaci")
lines.append(f"• Lista gena: preuzeta iz Cytoscape/STRING CSV ({os.path.basename(INPUT_CSV)}).")
lines.append(f"• Enrichment: g:Profiler (Enrichment_all u {os.path.basename(ENR_XLSX)}).")
lines.append(f"• Višestruko testiranje: adj_p_value (FDR), prag = {FDR_CUTOFF:.2f}.")
lines.append("")

lines.append("Rezultati (sažetak)")
lines.append(f"• Ukupan broj termina: {len(enr)}; značajnih na FDR ≤ {FDR_CUTOFF:.2f}: {n_sig}.")
lines.append(f"• Prikazani Top {min(TOP_N,len(top))} termina rangiranih po adj_p_value.")
lines.append("")

lines.append("Top termini")
for i, r in top.reset_index(drop=True).iterrows():
    glist = parse_intersections(r["intersections"])
    # dodaj smer uz gene, ako imamo logFC
    g_labels = []
    for g in glist:
        up = str(g).upper().strip()
        fc = logfc_map.get(up, None)
        if fc is None or pd.isna(fc):
            g_labels.append(g)
        else:
            arrow = "↑" if fc>0 else "↓" if fc<0 else "0"
            g_labels.append(f"{g} ({arrow} {fc:.2f})")
    line = (f"{i+1}. {r.get('Term','NA')} [{r.get('Source','NA')}] — "
            f"FDR={r['adj_p_value']:.2e}; hits={len(glist)}: {', '.join(g_labels)}")
    lines.append(line)

lines.append("")
lines.append("Napomena")
lines.append("• Ovo je automatski opis: ne sadrži tumačenja; navodi isključivo nazive termina, FDR i gene (sa znakom logFC).")
lines.append("• Za narativ po putu koristite Reactome/Metascape izveštaje; za AI sažetke uz kurirane izvore – QIAGEN IPA Interpret.")

os.makedirs(os.path.dirname(OUT_TXT), exist_ok=True)
with open(OUT_TXT, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print(f"[OK] Sačuvano: {OUT_TXT}")
