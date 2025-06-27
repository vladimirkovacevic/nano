from Bio import Entrez
from collections import defaultdict
from tqdm import tqdm
import time
import csv

# 🔧 Set your email for the Entrez API
Entrez.email = "your.email@example.com"

# 📚 Define search queries and corresponding output filenames
queries = [
    {
        "name": "microfluidics_pdms",
        "query": '''(
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
    },
    {
        "name": "polystyrene_nanoplastics",
        "query": '''(
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
    },
    {
        "name": "scrnaseq_pbmc",
        "query": '''
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
        "PBMC"[Title/Abstract] OR 
        "PBMCs"[Title/Abstract] OR 
        "peripheral blood mononuclear cells"[Title/Abstract] OR 
        "peripheral blood mononuclear cell"[Title/Abstract] OR 
        "human PBMCs"[Title/Abstract]
        )
        AND (review[Publication Type] OR journal article[Publication Type])
        '''
    },
    {
        "name": "scrnaseq_microfluidics",
        "query": '''(
        ("single cell RNA sequencing"[Title/Abstract] OR 
        "scRNA-seq"[Title/Abstract] OR 
        "single-cell RNA-seq"[Title/Abstract] OR 
        "single cell RNA-seq"[Title/Abstract] OR 
        "single cell transcriptomics"[Title/Abstract] OR 
        "single-cell transcriptomics"[Title/Abstract] OR 
        "single cell sequencing"[Title/Abstract])
        AND
        (
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
          OR
          "microfluidic chip"[Title/Abstract] OR 
          "microfluidic device"[Title/Abstract] OR 
          "lab-on-a-chip"[Title/Abstract] OR 
          "organ-on-a-chip"[Title/Abstract] OR 
          "microchannel"[Title/Abstract] OR 
          "single-cell microfluidics"[Title/Abstract] OR 
          "microfluidic platform"[Title/Abstract]
        )
        AND (review[Publication Type] OR journal article[Publication Type])
        )'''
    },
    {
        "name": "scrnaseq_title_only",
        "query": '''("single cell RNA sequencing"[Title] OR 
        "scRNA-seq"[Title] OR 
        "single-cell RNA-seq"[Title] OR 
        "single cell RNA-seq"[Title] OR 
        "single cell transcriptomics"[Title] OR 
        "single-cell transcriptomics"[Title] OR 
        "single cell sequencing"[Title])
        AND (review[Publication Type] OR journal article[Publication Type])'''
    },
    {
        "name": "scrnaseq_psnps",
        "query": '''
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
    },
]

# 📆 Define range of publication years
years = list(range(2001, 2026))

# 🔁 Loop through each query
for query_dict in queries:
    name = query_dict["name"]
    search_query = query_dict["query"]
    year_to_pmids = defaultdict(set)

    print(f"\n🔎 Running query: {name}")

    # 🔄 Loop through the years
    for year in tqdm(years, desc=f"Processing {name}"):
        date_query = f'("{year}"[Date - Publication])'
        final_query = f'{search_query} AND {date_query}'

        handle = Entrez.esearch(db="pubmed", term=final_query, retmax=100000)
        record = Entrez.read(handle)
        handle.close()

        pmids = set(record['IdList'])
        year_to_pmids[year] = pmids

        time.sleep(0.5)  # Be gentle to NCBI

    # 📊 Count unique PMIDs by year
    all_pmids = set()
    year_counts = []

    for year in years:
        unique_pmids = year_to_pmids[year] - all_pmids
        all_pmids.update(unique_pmids)
        year_counts.append((year, len(unique_pmids)))

    # 💾 Save to CSV file
    filename = f"{name}_filtered.csv"
    with open(filename, mode="w", newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Year", "Unique_PMID"])
        writer.writerows(year_counts)

    print(f"✅ Finished: Results saved to '{filename}'")

print("\n🎉 All queries processed successfully!")

