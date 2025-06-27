from Bio import Entrez
from collections import defaultdict
from tqdm import tqdm
import time
import csv

# 🔧 Postavi svoj email za Entrez API
Entrez.email = "tvoj.email@primer.com"

# 🔍 Upit za PBMCs
search_query = '''(
"PBMC"[Title/Abstract] OR 
"PBMCs"[Title/Abstract] OR 
"peripheral blood mononuclear cells"[Title/Abstract] OR 
"peripheral blood mononuclear cell"[Title/Abstract] OR 
"human PBMCs"[Title/Abstract]
)
AND (review[Publication Type] OR journal article[Publication Type])
'''

# 📆 Godine
years = list(range(2001, 2026))
year_to_pmids = defaultdict(set)

# 🔄 Petlja kroz godine
for year in tqdm(years, desc="Pretraga po godinama"):
    date_query = f'("{year}"[Date - Publication])'
    final_query = f'{search_query} AND {date_query}'

    handle = Entrez.esearch(db="pubmed", term=final_query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()

    pmids = set(record['IdList'])
    year_to_pmids[year] = pmids

    time.sleep(0.5)

# 📊 Broj unikatnih radova
all_pmids = set()
year_counts = []

for year in years:
    unique_pmids = year_to_pmids[year] - all_pmids
    all_pmids.update(unique_pmids)
    year_counts.append((year, len(unique_pmids)))

# 💾 Snimi CSV
with open("pbmcs_filtered.csv", mode="w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Godina", "Unikatni_PMID"])
    writer.writerows(year_counts)

print("✅ Završeno! Podaci su snimljeni u 'pbmcs_filtered.csv'.")
