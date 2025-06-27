from Bio import Entrez
from collections import defaultdict
from tqdm import tqdm
import time
import csv

# 🔧 OBAVEZNO: upiši svoj email
Entrez.email = "tvoj.email@primer.com"

# 🔍 Upit: scRNA-seq + polystyrene nanoplastics (samo naslovi i apstrakti)
search_query = '''
(
"single cell RNA sequencing"[Title/Abstract] OR 
"scRNA-seq"[Title/Abstract] OR 
"single-cell RNA-seq"[Title/Abstract] OR 
"single cell RNA-seq"[Title/Abstract] OR 
"single cell transcriptomics"[Title/Abstract] OR 
"single-cell transcriptomics"[Title/Abstract] OR 
"single cell sequencing"[Title/Abstract]
)
AND
(
"polystyrene nanoplastics"[Title/Abstract] OR 
"PSNPs"[Title/Abstract] OR 
"polystyrene nanoparticles"[Title/Abstract] OR 
"nano-polystyrene"[Title/Abstract] OR 
"polystyrene-based nanoparticles"[Title/Abstract] OR 
"polystyrene latex nanoparticles"[Title/Abstract] OR 
"PS nanoparticles"[Title/Abstract]
)
AND (review[Publication Type] OR journal article[Publication Type])
'''

# 📅 Godine
years = list(range(2001, 2026))
year_to_pmids = defaultdict(set)

# ⏳ Pretraga po godinama
for year in tqdm(years, desc="Preuzimanje po godinama"):
    date_query = f'("{year}"[Date - Publication])'
    final_query = f'{search_query} AND {date_query}'

    handle = Entrez.esearch(db="pubmed", term=final_query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()

    pmids = set(record['IdList'])
    year_to_pmids[year] = pmids
    time.sleep(0.5)

# 📊 Broj unikatnih PMID-a
all_pmids = set()
year_counts = []

for year in years:
    unique_pmids = year_to_pmids[year] - all_pmids
    all_pmids.update(unique_pmids)
    year_counts.append((year, len(unique_pmids)))

# 💾 Snimi CSV
with open("scRNAseq_PSNPs.csv", mode="w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Godina", "Unikatni_PMID"])
    writer.writerows(year_counts)

print("✅ Završeno! Fajl je snimljen kao scRNAseq_PSNPs.csv")
