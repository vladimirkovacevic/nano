from Bio import Entrez
from collections import defaultdict
from tqdm import tqdm
import time
import csv

# 🔧 Unesi svoj email za NCBI API
Entrez.email = "tvoj.email@primer.com"

# 🔍 Upit prema tvojoj modifikaciji
search_query = '''(
(
  "PDMS"[Title/Abstract] OR 
  "polydimethylsiloxane"[Title/Abstract] OR 
  "chip"[Title/Abstract]
)
AND 
(
  "microfluidic"[Title/Abstract] OR 
  "microfluidics"[Title/Abstract]
)
)
OR 
"microfluidic chip"[Title/Abstract] OR 
"microfluidic device"[Title/Abstract] OR 
"lab-on-a-chip"[Title/Abstract] OR 
"organ-on-a-chip"[Title/Abstract] OR 
"microchannel"[Title/Abstract] OR 
"single-cell microfluidics"[Title/Abstract] OR 
"microfluidic platform"[Title/Abstract]
AND (review[Publication Type] OR journal article[Publication Type])
'''

# 📅 Godine za pretragu
years = list(range(2001, 2026))
year_to_pmids = defaultdict(set)

# ⏳ Petlja po godinama
for year in tqdm(years, desc="Pretraga po godinama"):
    date_query = f'("{year}"[Date - Publication])'
    final_query = f'{search_query} AND {date_query}'

    handle = Entrez.esearch(db="pubmed", term=final_query, retmax=100000)
    record = Entrez.read(handle)
    handle.close()

    pmids = set(record['IdList'])
    year_to_pmids[year] = pmids

    time.sleep(0.5)  # poštuj NCBI pravila

# 📊 Brojanje unikatnih PMIDs
all_pmids = set()
year_counts = []

for year in years:
    unique_pmids = year_to_pmids[year] - all_pmids
    all_pmids.update(unique_pmids)
    year_counts.append((year, len(unique_pmids)))

# 💾 Snimi kao CSV
with open("microfluidics_PDMS_filtered.csv", mode="w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Godina", "Unikatni_PMID"])
    writer.writerows(year_counts)

print("✅ Završeno! Podaci su snimljeni u 'microfluidics_PDMS_filtered.csv'.")
