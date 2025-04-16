# SeuratX: R Package for zUMIs (Smart-seq3) Data Analysis and downstream analysis compatible with other scRNASeq techniques

**SeuratX** is an R package designed for streamlined quality control (QC), visualization, and analysis of Smart-seq3 single-cell RNA-seq data processed with **zUMIs**. It builds on top of the **Seurat** framework and introduces functions for reading expression matrices, computing QC metrics, and generating publication-ready plots.

---

## Installation

### 1. Install from GitHub

Make sure you have `devtools` installed:

```r
install.packages("devtools")
devtools::install_github("nimarafati/SeuratX")
```

> Replace `your-username` with your GitHub handle.

---

## Dependencies

- Seurat (>= 4.0)
- ggplot2
- gridExtra
- Matrix
- scales

Install any missing packages using:

```r
install.packages(c("Seurat", "ggplot2", "gridExtra", "Matrix", "scales"))
```

---

## Usage

### Step 1: Load metadata and read zUMIs outputs

`metadata` can look like this in `csv` format:
| plate_number | wells       | time   | part             | transgenic_line              | phenotype            | Genotype | Path                                                           | Sorting |
|--------------|-------------|--------|------------------|-------------------------------|----------------------|----------|----------------------------------------------------------------|---------|
| ZB_01_101    | full_plate  | 48_hpf | head_and_trunk   | kdrl:EGFP;gata1:dsRed         | enhanced_vascularity | WT1      | ~/projects/vascular_fate_study/data/ZB_01_101/                | sortA   |
| ZB_01_102    | whole_plate | 72_hpf | trunk_plus_tail  | fli1a:nls-mCherry;lyve1b:EGFP | lymphatic_shift      | KO1      | ~/experiments/lymph_shift/data/ZB_01_102/                     | sortB   |
| ZB_01_103    | quadrant_1  | 24_hpf | tail_bud         | lyve1b:GFP;mpeg1:mCherry      | no_phenotype         | MUT2     | ~/data/zebrafish_phenotyping/round3/ZB_01_103/                | sortC   |
| ZB_01_104    | full_plate  | 96_hpf | whole_body       | prox1a:GFP;kdrl:mCherry       | overgrowth_veins     | OE2      | ~/dev/zfish_screening_project/dataset_batch4/ZB_01_104/       | sortD   |
| ZB_01_105    | quadrant_2  | 30_hpf | somite_region    | flt4:YFP;sox17:dsRed          | angiogenic_delay     | WT2      | ~/labs/genomics_pipeline/test_runs/ZB_01_105/                 | sortE   |


```r
metadata <- 'smart-seq3_meta_data_MGKK_unique.txt'
proj_dir <- "your_project_path"

data_list <- SeuratX::read_zUMIs_data(
  metadata_file_name = metadata,
  metadata_dir = file.path(proj_dir, 'doc'),
  sample_name_column = 'plate_number',
  data_path_column = 'Path'
)
```

### Step 2: Merge data into a single Seurat object

```r
gtf_file <- "your_gtf_path.gtf"

res_list <- SeuratX::create_merged_seurat(
  data_list = data_list,
  feature_type = 'exon',
  gtf_file = gtf_file,
  gene_info_file = '../../data/zUMIs_annot.csv'
)

adata <- res_list$adata
genes.table <- res_list$genes.table
```

### Step 3: Compute QC statistics

```r
stats_list <- SeuratX::calculate_stats(
  adata = adata,
  filter_ercc = TRUE,
  filter_rRNA = TRUE,
  filter_mito = TRUE,
  mito_chr_name = "MT",
  genes_table = genes.table
)

adata <- stats_list$adata
genes.table <- stats_list$genes.table
```

### Step 4: Generate QC plots

```r
features <- c("nFeature_RNA", "percentMT", "percentERCC", "percentPC")
qc_dir <- file.path(proj_dir, 'results/QC')

vil_plt_list <- SeuratX::qc_plot(
  seurat_obj = adata,
  features = features,
  output_path = qc_dir,
  x_var = 'Plate',
  fill_var = 'Plate',
  show_cutoff = FALSE
)
```

### Step 5: Visualize relative expression

```r
SeuratX::compute_relative_expression(
  seurat_obj = adata,
  output_path = file.path(proj_dir, 'results/QC/'),
  top_gene = 20
)
```

### Step 6: Plot QC metrics on plate layout

```r
SeuratX::plot_qc_on_plate(
  adata = adata,
  features = features,
  output_dir = qc_dir,
  facet_var = 'Plate'
)
```

### Step 7: Filter cells by plate

```r
filt <- list()
for (plate in levels(as.factor(adata$orig.ident))) {
  filt[[plate]] <- SeuratX::cut_cells(adata@meta.data, plate = plate)
}
```

---

## Function Highlights

- `read_zUMIs_data()` — Load expression, UMI, gene counts, and barcodes.
- `create_merged_seurat()` — Combine all samples into a unified Seurat object.
- `calculate_stats()` — Compute ERCC, mitochondrial, rRNA percentages.
- `qc_plot()` — Generate violin plots with optional cutoffs.
- `compute_relative_expression()` — Profile top-expressed genes.
- `plot_qc_on_plate()` — Visualize features per well per plate.
- `cut_cells()` — Filter cells using thresholds or dynamic MAD filtering.

---

## 📁 Directory Structure Suggestion

```
project/
├── doc/
│   └── smart-seq3_meta_data_MGKK_unique.txt
├── data/
│   └── zUMIs_annot.csv
├── results/
│   └── QC/
├── scripts/
│   └── analysis.R
└── your_R_project.Rproj
```

---

## 📜 License

This package is released under the MIT License.

---

## 🤝 Contributions

Feel free to open issues, request features, or contribute improvements!
