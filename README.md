# Unraveling the molecular effects of nanoplastics on human blood cells using microfluidic chip-based single-cell RNA sequencing
To execute all analysis run `notebooks/full_analysis.ipynb`

## Dependencies

The `full_analysis.ipynb` notebook requires the following Python packages:

- `scanpy`
- `matplotlib`
- `pandas`
- `anndata`
- `numpy`
- `scipy`
- `seaborn`
- `gseapy`

You can install them using pip:

```bash
pip install scanpy matplotlib pandas anndata numpy scipy seaborn gseapy
```
## Data
scRNA reference cell types for PBMC (`cells.umi.new.txt, counts.umi.txt.gz, genes.umi.txt` and `meta.txt`) are downloaded from [Broad Institute repository](https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data#study-visualize) and merged into [sc.h5ad](https://drive.google.com/file/d/1u3bSTjpIuOwvHsOopd0tmtlO775dK-LV/view?usp=sharing)

Instructions for downloading sequenced samples are available in the analysis.
