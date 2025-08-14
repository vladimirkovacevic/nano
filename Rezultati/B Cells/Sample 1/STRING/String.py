# string_bcell_network.py
# Query STRING API za tvoje B-cell gene i eksport mreže za Cytoscape + brzi PNG pregled
# pip install requests pandas networkx matplotlib

import os, urllib.parse, requests, pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# === ULAZ/IZLAZ ===
CSV_IN   = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\B_cell_sample_1_ImmuneMeta_overlap.csv"
GENE_COL = "gene_name"
OUT_DIR  = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\STRING"

# === STRING PODESAVANJA ===
SPECIES     = 9606  # Homo sapiens
SCORE_CUTOFF= 0.4   # 0.15 (low) .. 0.4 (med) .. 0.7 (high)
STRING_BASE = "https://string-db.org/api"

def unique_genes_from_csv(path, col):
    df = pd.read_csv(path) if path.lower().endswith(".csv") else pd.read_excel(path)
    if col not in df.columns:
        raise SystemExit(f"[ERR] Ne postoji kolona '{col}' u {path}. Kolone: {list(df.columns)}")
    genes = (df[col].astype(str).str.strip().str.upper()
             .dropna().replace({"NAN": pd.NA}).dropna().unique().tolist())
    return genes

def string_get_ids(genes, species=9606):
    """Mapiranje simbola na STRING ID (/get_string_ids)"""
    ids_url = f"{STRING_BASE}/tsv/get_string_ids"
    identifiers = "%0d".join([urllib.parse.quote_plus(g) for g in genes])
    params = {"identifiers": identifiers, "species": species, "limit": 1}
    r = requests.get(ids_url, params=params, timeout=60)
    r.raise_for_status()
    rows, header = [], None
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if header is None:
            header = parts
            continue
        rows.append(dict(zip(header, parts)))
    return pd.DataFrame(rows)

def string_network(string_ids, species=9606, score_cutoff=0.4):
    """Preuzmi ivice mreže za zadate STRING ID-jeve (/network)"""
    net_url = f"{STRING_BASE}/tsv/network"
    identifiers = "%0d".join([urllib.parse.quote_plus(sid) for sid in string_ids])
    params = {"identifiers": identifiers, "species": species, "required_score": int(score_cutoff*1000)}
    r = requests.get(net_url, params=params, timeout=120)
    r.raise_for_status()
    rows, header = [], None
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if header is None:
            header = parts
            continue
        rows.append(dict(zip(header, parts)))
    return pd.DataFrame(rows)

def quick_preview(edges_df, out_png):
    """Brzi networkx preview (samo ilustracija; za publikaciju koristi Cytoscape)"""
    if edges_df.empty:
        print("[WARN] Nema ivica za preview."); return
    G = nx.Graph()
    ucol = "preferredName_A" if "preferredName_A" in edges_df else "stringId_A"
    vcol = "preferredName_B" if "preferredName_B" in edges_df else "stringId_B"
    wcol = "score" if "score" in edges_df else None

    for _, row in edges_df.iterrows():
        u, v = row[ucol], row[vcol]
        w = float(row[wcol]) if wcol and str(row[wcol]).replace(".","",1).isdigit() else 0.0
        G.add_edge(u, v, weight=w)

    pos = nx.spring_layout(G, seed=42, k=0.6)
    plt.figure(figsize=(8, 6))
    nx.draw_networkx_nodes(G, pos, node_size=600)
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[OK] PNG preview: {out_png}")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    genes = unique_genes_from_csv(CSV_IN, GENE_COL)
    print(f"[INFO] Učitano gena: {len(genes)}")
    if not genes:
        raise SystemExit("[ERR] Lista gena je prazna.")

        # 1) Mapiranje na STRING ID
    map_df = string_get_ids(genes, SPECIES)
    if map_df.empty:
        raise SystemExit("[ERR] STRING mapping vratio prazan skup.")

    # obezbedi da postoji kolona 'score'
    if "score" in map_df.columns:
        map_df["score"] = pd.to_numeric(map_df["score"], errors="coerce").fillna(0.0)
    else:
        map_df["score"] = 0.0

    # preimenuj kolone pre sortiranja
    if "queryItem" in map_df.columns:
        map_df.rename(columns={"queryItem": "input_gene"}, inplace=True)
    if "preferredName" in map_df.columns:
        map_df.rename(columns={"preferredName": "preferredName"}, inplace=True)  # no-op, ali ostavimo radi jasnoće

    # zadrži best match po input genu
    sort_cols = ["input_gene", "score"] if "input_gene" in map_df.columns else ["score"]
    map_df = map_df.sort_values(sort_cols, ascending=[True, False] if len(sort_cols) == 2 else False)
    map_df = map_df.drop_duplicates("input_gene" if "input_gene" in map_df.columns else "stringId")


    # 2) Mreža
    edges_df = string_network(map_df["stringId"].tolist(), SPECIES, SCORE_CUTOFF)

    # 3) Dodaj nice names
    id2name = dict(zip(map_df["stringId"], map_df["preferredName"]))
    edges_df["preferredName_A"] = edges_df["stringId_A"].map(id2name)
    edges_df["preferredName_B"] = edges_df["stringId_B"].map(id2name)

    # 4) Sačuvaj fajlove
    edges_out = os.path.join(OUT_DIR, "B_cell_sample_1_STRING_edges.tsv")
    nodes_out = os.path.join(OUT_DIR, "B_cell_sample_1_STRING_nodes.tsv")
    map_df.to_csv(nodes_out, sep="\t", index=False)
    edges_df.to_csv(edges_out, sep="\t", index=False)
    print(f"[OK] Sačuvano:\n  {nodes_out}\n  {edges_out}")

    # 5) PNG preview
    png_out = os.path.join(OUT_DIR, "B_cell_sample_1_STRING_preview.png")
    quick_preview(edges_df, png_out)

if __name__ == "__main__":
    main()
