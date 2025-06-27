from Bio import Entrez
from collections import defaultdict
from tqdm import tqdm
import time
import csv

# 🔧 Podesi svoj email
Entrez.email = "tvoj.email@primer.com"

# 🔍 Pretraga - isključivo u naslovu
search_query = '''("single cell RNA sequencing"[Title] OR 
"scRNA-seq"[Title] OR 
"single-cell RNA-seq"[Title] OR 
"single cell RNA-seq"[Title] OR 
"single cell transcriptomics"[Title] OR 
"single-cell transcriptomics"[Title] OR 
"single cell sequencing"[Title])
AND (review[Publication Type] OR journal article[Publication Type])'''

# 📆 Godine
years = list(range(2001, 2026))
year_to_pmids = defaultdict(set)

# 🔄 Pretraga po godinama
for year in tqdm(years, desc="Preuzimanje po godinama"):
    date_query = f'("{year}"[Date - Publication])'
    final_query = f'{search_query} AND {date_query}'

    handle = Entrez.esearch(db="pubmed", term=final_query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()

    pmids = set(record['IdList'])
    year_to_pmids[year] = pmids

    time.sleep(0.5)

# 📊 Brojanje unikatnih radova
all_pmids = set()
year_counts = []

for year in years:
    unique_pmids = year_to_pmids[year] - all_pmids
    all_pmids.update(unique_pmids)
    year_counts.append((year, len(unique_pmids)))

# 💾 Snimi CSV
with open("scRNAseq_title_only.csv", mode="w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Godina", "Unikatni_PMID"])
    writer.writerows(year_counts)

print("✅ Završeno! Podaci su snimljeni u 'scRNAseq_title_only.csv'.")
