import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# ================== PUTANJE (prilagođeno tvom fajlu) ==================
INPUT_CSV = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Cytoscape\B_cell_sample_1_ImmuneMeta_overlap_STRING.csv"
OUT_DIR   = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Bar Plot Up-DownRegulation"
OUT_PNG   = os.path.join(OUT_DIR, "B_cell_sample_1_ImmuneMeta_overlap_BAR.png")

# Koliko gena crtati (None = sve; npr. 40 za top 40 po |logFC|)
TOP_N = None

# ================== POMOĆNE ==================
def sniff_sep(path):
    """Vrati najverovatniji separator iz prve 2–3 linije."""
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        head = "".join([next(f, "") for _ in range(3)])
    cands = [",", "\t", ";", "|"]
    counts = {sep: head.count(sep) for sep in cands}
    sep = max(counts, key=counts.get)
    return sep if counts[sep] > 0 else ","

def find_col(candidates, cols):
    """Nađi kolonu neosetljivo na velika/mala slova i razmake."""
    norm = {re.sub(r"\s+", "", c).lower(): c for c in cols}
    for c in candidates:
        key = re.sub(r"\s+", "", c).lower()
        if key in norm:
            return norm[key]
    return None

# ================== UČITAVANJE ==================
if not os.path.exists(INPUT_CSV):
    raise SystemExit(f"[GREŠKA] Ne nalazim ulazni fajl:\n{INPUT_CSV}")

sep = sniff_sep(INPUT_CSV)
df = pd.read_csv(INPUT_CSV, sep=sep, encoding="utf-8-sig")

# ================== KOLONE ==================
gene_col  = find_col(["gene_name", "gene", "symbol"], df.columns)
logfc_col = find_col(["logfoldchanges", "log2fc", "logfc"], df.columns)

if gene_col is None or logfc_col is None:
    raise SystemExit(f"[GREŠKA] Ne nalazim obavezne kolone (gene_name / logfoldchanges)."
                     f"\nUčitane kolone su: {list(df.columns)}\nSeparator detektovan: {repr(sep)}")

# ================== PRIPREMA ==================
# Sortiranje (po potrebi top N po apsolutnoj promeni)
if TOP_N:
    df = df.reindex(df[logfc_col].abs().sort_values(ascending=False).index).head(TOP_N)
else:
    df = df.sort_values(logfc_col, ascending=False)

colors = ["red" if v > 0 else "blue" for v in df[logfc_col]]

# Dinamička širina figure da se vide labele
n = len(df)
fig_w = max(10, min(0.35 * n, 50))

# ================== PLOT ==================
plt.figure(figsize=(fig_w, 6))
plt.bar(df[gene_col].astype(str), df[logfc_col], color=colors)
plt.axhline(0, color="black", linewidth=0.8)
plt.xlabel("Gene")
plt.ylabel(logfc_col)
plt.title("Up- i down-regulisani geni (B cell sample 1 – Immune Meta)")
plt.xticks(rotation=75, ha="right")
plt.tight_layout()

os.makedirs(OUT_DIR, exist_ok=True)
plt.savefig(OUT_PNG, dpi=600)
plt.show()

print(f"[INFO] Sačuvano: {OUT_PNG}\n[INFO] Korišćen separator: {repr(sep)}")
