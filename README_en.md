<!-- Language Toggle: English README，按键切回中文主 README -->
<div id="lang-toggle" style="margin-bottom: 12px;">
  <a href="./README.md">
    <button>中文</button>
  </a>
  <a href="./README_en.md" style="margin-left: 8px;">
    <button>English</button>
  </a>
</div>

> This is the **English README**. Click the **中文** button above to go back to the main Chinese README.

## Project Overview: Pheno Fusion GenNet Multitask Model

PyTorch-based multitask interaction model with optional GenNet topology (SNP→Gene) for **genotype × environment × multiple phenotypes**. Canonical entry: `main_gennet.py`.

---

## Environments
- Python = 3.9
- R = 4.4.2
- OS: Windows / Linux
- Python deps: see `requirements.txt`

---

## Project Layout (key items)
- `main_gennet.py` — main training/_eval entry
- `hyperparameter_search_v3.py` — optional hyperparam search
- `data/` — aligned inputs (`aligned_GENO.*`, `aligned_ECOV.*`, `South_PHENO.*`)
- `data_process/data_aligned.py` — CSV↔Parquet + alignment loader
- `models/` — core model `interaction_model_gennet.py`, encoders, GenNet layer
- `results/` — training plots; `visualization_results/` — sample visuals

---

## End-to-End Workflow (R → Python)

### Step 1. R setup
Install required R packages:
```r
install.packages(c(
  "tidyverse", "data.table", "lubridate", "jsonlite",
  "apsimx", "soilDB", "spData", "sf",
  "ggrepel"
))
```

### Step 2. Raw data curation (optional, from G2F)
In `data_curation_and_ecov/`, run modules in order:
```r
setwd("data_curation_and_ecov")
source("1_phenotypes/code/1_get_phenotype.R")          # -> PHENO.csv
# (HPC) run 1_SNP_filter_codify.sub, then:
source("2_genotypes/code/2_LD_prune_genotypes.R")      # -> GENO.csv, MAP.csv
source("3_year_loc_summary/code/1_year_location_summary.R")
source("4_weather/code/1_get_met_apsim.R")             # optional
source("5_APSIM/code/1_apsim_sim.R")                   # optional, heavy
source("6_ecov/code/1_get_env_cov.R")                  # -> ECOV.csv
source("7_GDD/code/1_add_GDD_to_pheno.R")              # add GDD traits
```
If you already have PHENO/GENO/ECOV, you can skip this step.

### Step 3. Data alignment (Python)
1) Create env and install deps:
```bash
conda create -n GxE_env python==3.12.0
conda activate GxE_env
pip install -r requirements.txt
```
2) Adjust paths in `data_process/data_aligned.py`, then run:
```bash
python data_process/data_aligned.py
```
3) Place aligned files into `data/`:
- `data/aligned_GENO.(csv|parquet)`
- `data/aligned_ECOV.(csv|parquet)`
- `data/South_PHENO.(csv|parquet)` (includes the 7 target traits)

### Step 4. Train
- Basic (no topology):
```bash
python main_gennet.py --epochs 100
```
- With GenNet topology:
```bash
python main_gennet.py \
  --use_gennet \
  --topology_file gennet_data/topology_gene.csv \
  --epochs 120
```
- With hyperparam search (optional):
```bash
python main_gennet.py --hyperparameter_search --epochs 100
```

### Step 5. Outputs
- Checkpoints: `final_model_gennet_YYYYMMDD_HHMMSS.pth`
- Plots in `results/`: training loss, per-trait metrics, task relationships/weights
- Extra visuals in `visualization_results/`

---

## Inputs expected by `main_gennet.py`
- Genotype: `data/aligned_GENO.(csv|parquet)` (first column dropped)
- Environment: `data/aligned_ECOV.(csv|parquet)` (first column dropped)
- Phenotype: `data/South_PHENO.(csv|parquet)` with 7 target columns
- All three tables must be row-aligned (same samples, same order).

---

## Run Quick Checklist
- R 4.4.2 + required R packages (if running the R pipeline)
- Python env active; `pip install -r requirements.txt`
- `data/` contains aligned GENO/ECOV/PHENO with matching rows
- (Optional) topology CSV ready if using `--use_gennet`


