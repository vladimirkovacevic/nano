import pandas as pd
import requests
import os

# Ulazni CSV sa listom gena (12 gena nakon filtriranja)
IN_CSV = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Pathway Enrichment\B_cell_sample_1_ImmuneMeta_overlap_STRING.csv"

# Kolona sa imenima gena
GENE_COL = "gene_name"

# Folder za rezultate
OUT_DIR = os.path.dirname(IN_CSV)
OUT_FILE = os.path.join(OUT_DIR, "B_cell_sample_1_PathwayEnrichment_Enrichr.csv")

# Proveri da li fajl postoji
if not os.path.exists(IN_CSV):
    raise SystemExit(f"[GREŠKA] Ulazni fajl ne postoji: {IN_CSV}")

# Učitaj listu gena
df = pd.read_csv(IN_CSV)
genes = df[GENE_COL].dropna().unique().tolist()

print(f"[INFO] Učitano gena: {len(genes)}")

# Enrichr API – set za analizu
ENRICHR_LIBRARIES = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "GO_Biological_Process_2023"
]

def enrichr_post_list(genes):
    url = "https://maayanlab.cloud/Enrichr/addList"
    payload = {
        "list": "\n".join(genes),
        "description": "B_cell_sample_1 pathway enrichment"
    }
    response = requests.post(url, files=payload)
    return response.json()

def enrichr_get_results(user_list_id, library):
    url = f"https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={library}"
    response = requests.get(url)
    return response.json()

# Pošalji listu gena
print("[INFO] Slanje liste gena na Enrichr...")
res = enrichr_post_list(genes)
user_list_id = res.get("userListId")
if not user_list_id:
    raise SystemExit("[GREŠKA] Nije vraćen userListId.")

# Skupi rezultate iz svih biblioteka
all_results = []
for lib in ENRICHR_LIBRARIES:
    print(f"[INFO] Analiza biblioteke: {lib}")
    results = enrichr_get_results(user_list_id, lib)
    for r in results.get(lib, []):
        term_name = r[1]
        p_value = r[2]
        combined_score = r[4]
        genes_in_term = r[5]
        all_results.append({
            "Library": lib,
            "Term": term_name,
            "P_value": p_value,
            "Combined_score": combined_score,
            "Genes": genes_in_term
        })

# Snimi rezultate
pd.DataFrame(all_results).to_csv(OUT_FILE, index=False)
print(f"[GOTOVO] Rezultati snimljeni: {OUT_FILE}")
