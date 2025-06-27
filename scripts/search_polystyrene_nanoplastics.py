from Bio import Entrez
from collections import defaultdict
from tqdm import tqdm
import time
import csv

# 🔧 Neophodno za NCBI API
Entrez.email = "tvoj.email@primer.com"  # ← stavi svoj email

# 🔍 Sinonimi za polystyrene nanoplastics
search_query = '''(
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

# 📅 Godine koje analiziramo
years = list(range(2001, 2026))
year_to_pmids = defaultdict(set)

# ⏳ Iteracija kroz godine
for year in tqdm(years, desc="Pretraga po godinama"):
    date_query = f'("{year}"[Date - Publication])'
    final_query = f'{search_query} AND {date_query}'

    handle = Entrez.esearch(db="pubmed", term=final_query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()

    pmids = set(record['IdList'])
    year_to_pmids[year] = pmids

    time.sleep(0.5)  # poštuj NCBI pravila

# 📊 Broj unikatnih radova po godini
all_pmids = set()
year_counts = []

for year in years:
    unique_pmids = year_to_pmids[year] - all_pmids
    all_pmids.update(unique_pmids)
    year_counts.append((year, len(unique_pmids)))

# 💾 Snimanje u CSV
with open("psnps_filtered.csv", mode="w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Godina", "Unikatni_PMID"])
    writer.writerows(year_counts)

print("✅ Završeno! Podaci su snimljeni u 'psnps_filtered.csv'.")
