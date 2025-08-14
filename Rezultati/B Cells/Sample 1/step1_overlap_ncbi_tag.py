# -*- coding: utf-8 -*-
# step1_overlap_ncbi_tag.py
#
# Radi 3 stvari:
#  1) Učita DE fajl (B_cell_sample_1.csv) i 6 Hallmark .gmt setova → napravi presek (overlap)
#  2) Izračuna prebrojavanja po setovima
#  3) Za overlapp gene povuče NCBI summary (E-utilities) i dodeli "molecule_class"
#
# Izlazi (u BASE_DIR):
#  - B_cell_sample_1_ImmuneMeta_overlap.csv
#  - B_cell_sample_1_ImmuneMeta_counts_by_set.csv
#  - B_cell_sample_1_Tagged_All.csv

import os, time
from collections import defaultdict
import pandas as pd
import requests

# =====================
# PODEŠAVANJA PUTANJA
# =====================
BASE_DIR = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1"
# Umesto stare lokacije, direktno postavljamo apsolutnu putanju na HALLMARKS folder:
HALLMARKS_DIR = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\Sources\HALLMARKS"
DE_PATH = os.path.join(BASE_DIR, "B_cell_sample_1.csv")

# .gmt fajlovi koje želiš da koristiš
HALLMARK_FILES = [
    "HALLMARK_INFLAMMATORY_RESPONSE.v2025.1.Hs.gmt",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB.v2025.1.Hs.gmt",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING.v2025.1.Hs.gmt",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE.v2025.1.Hs.gmt",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE.v2025.1.Hs.gmt",
    "HALLMARK_COMPLEMENT.v2025.1.Hs.gmt",
]
HALLMARK_PATHS = [os.path.join(HALLMARKS_DIR, f) for f in HALLMARK_FILES]

# Izlazni fajlovi
OUT_OVERLAP = os.path.join(BASE_DIR, "B_cell_sample_1_ImmuneMeta_overlap.csv")
OUT_COUNTS  = os.path.join(BASE_DIR, "B_cell_sample_1_ImmuneMeta_counts_by_set.csv")
OUT_TAGGED  = os.path.join(BASE_DIR, "B_cell_sample_1_Tagged_All.csv")

# NCBI podešavanje
NCBI_EMAIL = "marko.zivanovic@uni.kg.ac.rs"
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# =====================
# POMOĆNE FUNKCIJE
# =====================
def find_col(cands, cols):
    low = [c.lower() for c in cols]
    for c in cands:
        if c.lower() in low:
            return cols[low.index(c.lower())]
    return None

def load_meta_set(gmt_paths):
    """Učitava više .gmt fajlova i vraća DataFrame: gene, in_sets"""
    for p in gmt_paths:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Nedostaje .gmt fajl: {p}")

    gene_to_sets = defaultdict(set)
    present_names = set()
    for gmt in gmt_paths:
        with open(gmt, "r", encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    parts = line.split()  # fallback ako su razmaci
                if len(parts) < 3:
                    continue
                set_name = parts[0].strip()
                present_names.add(set_name)
                genes = [g.strip().upper() for g in parts[2:] if g.strip()]
                for g in genes:
                    gene_to_sets[g].add(set_name)

    meta = pd.DataFrame(
        [{"gene": g, "in_sets": ";".join(sorted(s))} for g, s in gene_to_sets.items()]
    ).sort_values("gene").reset_index(drop=True)
    return meta, sorted(present_names)

def ncbi_gene_summary(gene, email=NCBI_EMAIL):
    """Vrati NCBI summary za dati gen (Homo sapiens)."""
    from urllib.parse import quote_plus
    try:
        term = f'{gene}[Gene Name] AND "Homo sapiens"[Organism]'
        esearch = f"{NCBI_BASE}/esearch.fcgi?db=gene&term={quote_plus(term)}&retmode=json&email={email}"
        r = requests.get(esearch, timeout=30)
        if r.status_code != 200:
            return {}
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return {}
        gid = ids[0]
        esum = f"{NCBI_BASE}/esummary.fcgi?db=gene&id={gid}&retmode=json&email={email}"
        r2 = requests.get(esum, timeout=30)
        if r2.status_code != 200:
            return {"ncbi_geneid": gid}
        res = r2.json().get("result", {}).get(gid, {})
        out = {
            "ncbi_geneid": gid,
            "ncbi_name": res.get("name"),
            "ncbi_description": res.get("description"),
            "ncbi_summary": res.get("summary")
        }
        if out.get("ncbi_summary"):
            out["ncbi_summary_short"] = out["ncbi_summary"].replace("\n", " ")[:300]
        return out
    except Exception:
        return {}

# Heurističko tagovanje
CANON_TF = {"NFKB1","NFKB2","RELA","RELB","REL","JUN","FOS","STAT1","STAT3","STAT5A","STAT5B","IRF1","IRF7","IRF3"}
CYTOKINES = {"IL1B","IL6","IL8","TNF","IFNG","IFNA1","IFNA2"}
CHEMOKINE_PREFIXES = ("CXCL","CCL","CX3CL","XCL")
RECEPTOR_SUFFIXES  = ("R","RA","RB")

def tag_class(symbol, summary_text=""):
    s = str(symbol).strip().upper()
    t = (summary_text or "").lower()
    def has(*kw): return any(k.lower() in t for k in kw)

    if s in CANON_TF or "transcription factor" in t or "dna-binding transcription factor" in t:
        return "transcription factor"
    if s in CYTOKINES or s.startswith(CHEMOKINE_PREFIXES) or has("cytokine","chemokine"):
        return "signaling molecule (cytokine/chemokine)"
    if s.endswith(RECEPTOR_SUFFIXES) or has("receptor"):
        return "receptor"
    if has("kinase","jak","ikk","tbk1","map kinase"):
        return "enzyme (kinase)"
    if has("toll-like receptor","tlr","pattern recognition","rig-i","cgas","sting","nod-like"):
        return "pattern-recognition receptor / innate sensor"
    if s == "NFKBIA" or has("inhibitor","negative regulator","ikb","iκb"):
        return "inhibitor / negative regulator"
    if has("adaptor","adapter","scaffold","traf","myd88"):
        return "adaptor / scaffold"
    return "other/immune"

# =====================
# GLAVNI TOK
# =====================
def main():
    os.makedirs(BASE_DIR, exist_ok=True)

    # 1) Učitaj DE
    de = pd.read_csv(DE_PATH)
    gene_col = find_col(["gene_name","gene","symbol"], de.columns)
    if not gene_col:
        raise SystemExit(f"[GREŠKA] Nema kolone gena u DE fajlu. Kolone su: {list(de.columns)}")
    de["_gene_upper"] = de[gene_col].astype(str).str.strip().str.upper()

    # 2) Meta-set iz .gmt fajlova
    meta_df, present_sets = load_meta_set(HALLMARK_PATHS)
    meta_df["_gene_upper"] = meta_df["gene"].astype(str).str.upper()

    # 3) Presek (overlap)
    overlap = de.merge(meta_df[["_gene_upper","in_sets"]], on="_gene_upper", how="inner")

    # 4) Statistika po setovima
    counts_total = {s: (meta_df["in_sets"].str.contains(s).sum()) for s in present_sets}
    counts_overlap = {s: (overlap["in_sets"].str.contains(s).sum()) for s in present_sets}
    counts_df = pd.DataFrame({
        "set": present_sets,
        "genes_in_set_total": [counts_total[s] for s in present_sets],
        "genes_in_overlap": [counts_overlap[s] for s in present_sets],
    })
    counts_df["overlap_fraction_%"] = (
        counts_df["genes_in_overlap"] / counts_df["genes_in_set_total"].replace(0, pd.NA)
    ) * 100

    # 5) NCBI anotacija za jedinstvene overlapp gene
    genes = overlap[gene_col].astype(str).str.upper().unique().tolist()
    ann_rows = []
    for g in genes:
        info = ncbi_gene_summary(g)
        info["gene"] = g
        ann_rows.append(info)
        time.sleep(0.34)   # NCBI rate-limit
    ann = pd.DataFrame(ann_rows)

    # 6) Merge + tagovanje
    overlap["_GENE_UPPER_"] = overlap[gene_col].astype(str).str.upper()
    if not ann.empty:
        ann["_GENE_UPPER_"] = ann["gene"].astype(str).str.upper()
        overlap = overlap.merge(
            ann[["_GENE_UPPER_","ncbi_summary","ncbi_summary_short"]],
            on="_GENE_UPPER_", how="left"
        )
    else:
        overlap["ncbi_summary"] = pd.NA
        overlap["ncbi_summary_short"] = pd.NA

    overlap["molecule_class"] = [
        tag_class(sym, summ) for sym, summ in zip(overlap[gene_col], overlap["ncbi_summary_short"])
    ]
    overlap.drop(columns=["_GENE_UPPER_"], inplace=True, errors="ignore")

    # 7) Snimi rezultate
    overlap.to_csv(OUT_OVERLAP, index=False, encoding="utf-8-sig")
    counts_df.to_csv(OUT_COUNTS, index=False, encoding="utf-8-sig")

    # Sažeta tabela za rad u Cytoscape itd.
    fc_col   = find_col(["logfoldchanges","log2fc","logfc"], overlap.columns)
    absfc_col= find_col(["logfoldchanges_abs","abs_log2fc","abs_logfc"], overlap.columns)
    padj_col = find_col(["pvals_adj","padj","qval","fdr","pval_adj"], overlap.columns)

    keep = [gene_col, fc_col, absfc_col, padj_col, "in_sets", "molecule_class", "ncbi_summary_short"]
    keep = [c for c in keep if c]  # izbaci None
    tagged = overlap[keep].copy()
    tagged = tagged.rename(columns={
        gene_col: "gene",
        fc_col: "log2FC",
        absfc_col: "|log2FC|",
        padj_col: "FDR",
        "ncbi_summary_short": "NCBI_summary_short",
    })
    tagged.to_csv(OUT_TAGGED, index=False, encoding="utf-8-sig")

    print("— GOTOVO —")
    print(f"DE zapisa: {len(de)} | meta-set gena: {len(meta_df)} | presek: {len(overlap)}")
    print("Snimljeno:")
    print(" ", OUT_OVERLAP)
    print(" ", OUT_COUNTS)
    print(" ", OUT_TAGGED)

if __name__ == "__main__":
    main()
