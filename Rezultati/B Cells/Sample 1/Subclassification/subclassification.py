import pandas as pd
import os

# === PATH CONFIG ===
base_dir = r"D:\1 Rad\BioINGLab\scRNA-Seq\Marko Genes Analyses\Rezultati\B Cells\Sample 1\Subclassification"
input_file = os.path.join(base_dir, "Gene description.txt")
output_file = os.path.join(base_dir, "B_cell_sample_1_ImmuneMeta_overlap_subclassed.csv")

# === Load Data ===
df = pd.read_csv(input_file, sep="\t", encoding="utf-8")  # ako nije tab, promeni na sep=","

# === Function for subclassification ===
def classify(summary):
    if pd.isna(summary):
        return "Unclassified"
    s = summary.lower()

    if "ubiquitin" in s and "ligase" in s:
        return "Ubiquitin ligase / protein degradation regulator"
    elif "interferon" in s or "antiviral" in s:
        return "Interferon-stimulated antiviral protein"
    elif "glycosyltransferase" in s or "glycosylation" in s:
        return "Glycosyltransferase / cell adhesion modulator"
    elif "phosphatase" in s:
        return "Signaling modulator (phosphatase)"
    elif "kynuren" in s or "tryptophan" in s:
        return "Tryptophan metabolism enzyme (kynurenine pathway)"
    elif "tnf" in s or "tumor necrosis" in s or "ligand" in s:
        return "TNF family member / inflammatory signal"
    elif "hla" in s or "major histocompatibility" in s:
        return "HLA antigen presentation"
    elif "chemokine" in s or "cytokine" in s:
        return "Cytokine / chemokine signaling"
    elif "receptor" in s:
        return "Immune receptor"
    else:
        return "Other immune-related function"

# === Apply only to "other/immune" ===
df["molecule_subclass"] = df.apply(
    lambda row: classify(row["ncbi_summary"]) if row["molecule_class"] == "other/immune" else row["molecule_class"],
    axis=1
)

# === Save output ===
df.to_csv(output_file, index=False, encoding="utf-8")

print(f"✅ Subclassification complete. Saved to:\n{output_file}")
