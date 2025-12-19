<!-- Language Toggle -->
<div id="lang-toggle" style="margin-bottom: 12px;">
  <button onclick="document.getElementById('lang-zh').style.display='block';document.getElementById('lang-en').style.display='none';">
    中文
  </button>
  <button onclick="document.getElementById('lang-en').style.display='block';document.getElementById('lang-zh').style.display='none';">
    English
  </button>
</div>

<div id="lang-zh" style="display: block;">

## 项目简介：Pheno Fusion GenNet 多任务表型预测模型

本项目基于 PyTorch 实现了一个 **集成 GenNet 拓扑结构的多任务交互模型**，用于玉米等作物的 **基因型（SNP）×环境（ECOV）×多性状表型** 联合建模与预测。  
当前代码已精简为以 `main_gennet.py` 为**唯一训练入口**，并保留了所有必要的模型与数据处理组件。

- **多任务学习**：同时预测多个农艺性状（如 `yield / silking / anthesis / grain_moisture` 等）
- **GenNet 集成**：可选使用 SNP→基因拓扑（基于 `topology_file`），提升可解释性
- **任务关系学习**：显式学习任务间相关性矩阵，并可视化任务关系与任务权重

---

## 环境依赖

- Python = 3.9
- R = 4.4.2
- 操作系统：Windows/Linux
- 主要依赖请参照 `requirements.txt`。

---

## 项目结构

- `main_gennet.py`：**项目主入口**，完成数据加载、划分、模型构建、训练、评估与可视化、模型保存
- `hyperparameter_search_v3.py`：多任务模型的超参数搜索脚本（可选，通过 `--hyperparameter_search` 启用）
- `train_eval.py`：早期版本的训练 / 评估函数（当前入口脚本未直接调用，但仍保留以便复用）
- `data/`：基础数据文件
  - `aligned_GENO.(csv|parquet)`：对齐后的 SNP 基因型矩阵
  - `aligned_ECOV.(csv|parquet)`：对齐后的环境协变量（ECOV）矩阵
  - `South_PHENO.(csv|parquet)`：南方试验点的表型数据（包含多个目标性状列）
  - 其他如 `ECOV_KL.csv`, `MAP.csv`, `PHENO.csv` 等为扩展/原始数据，可按需使用
- `data_process/data_aligned.py`：封装 `load_or_convert_data`，支持 CSV→Parquet 自动转换与加载
- `models/`：模型与子模块
  - `interaction_model_gennet.py`：**核心多任务 GenNet 模型 `MultiTaskGenNetModel`**
  - `interaction_model_v3.py` ：多任务交互模型
  - `encoders_v2.py`：基因型与环境编码器、多头自注意力模块
  - `feature_selector.py`：稀疏特征选择层（SNP gating）
  - `gennet_layer.py`：GenNet 的 `LocallyDirected1D` 稀疏拓扑层与拓扑 mask 构建工具
- `results/`：训练过程中自动保存的图像结果（损失曲线、任务关系热力图、任务权重柱状图等）
- `visualization_results/`：额外的任务关系可视化示例图

---

## 数据准备

### 0. 原始数据处理 Pipeline（R 语言）

在进入 Python 训练流程之前，需要从原始 G2F（Genomes to Fields）数据生成标准化的表型、基因型、环境协变量文件。本项目提供了完整的 R 语言数据处理 pipeline，位于 `data_curation_and_ecov/` 目录。

#### Pipeline 结构

该 pipeline 包含 **7 个顺序执行的模块**，每个模块的代码位于对应文件夹的 `code/` 子目录，输出文件保存在 `output/` 子目录：

```
data_curation_and_ecov/
├── 1_phenotypes/          # 模块1：表型数据整理
│   └── code/1_get_phenotype.R
├── 2_genotypes/          # 模块2：基因型数据过滤与编码
│   ├── code/1_SNP_filter_codify.sub
│   └── code/2_LD_prune_genotypes.R
├── 3_year_loc_summary/   # 模块3：年份-地点摘要
│   └── code/1_year_location_summary.R
├── 4_weather/            # 模块4：天气数据下载
│   └── code/1_get_met_apsim.R
├── 5_APSIM/              # 模块5：APSIM 作物模型模拟
│   └── code/1_apsim_sim.R
├── 6_ecov/               # 模块6：环境协变量生成
│   └── code/1_get_env_cov.R
├── 7_GDD/                # 模块7：GDD 性状计算
│   └── code/1_add_GDD_to_pheno.R
├── source/               # 原始数据（从 G2F 网站下载）
│   ├── Phenotype/
│   ├── Genotype/
│   ├── Metadata/
│   └── Agronomic_information/
└── tools/                # R 工具函数库
    ├── Functions.R
    ├── read_phenotype.R
    ├── read_metadata.R
    ├── LD_prune.R
    ├── ecov_utils.R
    └── APSIM_functions.R
```

#### 数据来源

原始数据需要从 **G2F 网站**（https://www.genomes2fields.org/resources/）下载，并按以下结构放置在 `source/` 目录：

- **`source/Phenotype/`**：各年份的表型数据 CSV 文件（如 `g2f_2014_hybrid_data_clean.csv`, `g2f_2015_hybrid_data_clean.csv`, ...）
- **`source/Genotype/`**：基因型 VCF 文件（如 `Hybrids_G2F_Genotype_Data_All_Years.vcf.zip`）
- **`source/Metadata/`**：各年份的试验地点元数据（如 `g2f_2014_field_characteristics.csv`, `g2f_2015_field_metadata.csv`, ...）
- **`source/Agronomic_information/`**：农艺信息文件（如 `g2f_2014_field_characteristics.csv`, `g2f_2015_agronomic_information.csv`, ...）

#### 各模块处理流程

**模块 1：表型数据整理** (`1_phenotypes/code/1_get_phenotype.R`)
- 读取各年份的表型数据、元数据和农艺信息
- 合并表型与元数据，统一基因型名称格式
- 输出：`1_phenotypes/output/PHENO.csv`

**模块 2：基因型数据处理** (`2_genotypes/`)
- **步骤 1** (`1_SNP_filter_codify.sub`)：从 VCF 文件提取 SNP，进行 MAF（最小等位基因频率）和缺失率过滤，编码为 0/1/2 格式
- **步骤 2** (`2_LD_prune_genotypes.R`)：对 SNP 进行 LD（连锁不平衡）剪枝，减少冗余位点
- 输出：`2_genotypes/output/GENO.csv`（基因型矩阵）和 `MAP.csv`（SNP 位置信息）

**模块 3：年份-地点摘要** (`3_year_loc_summary/code/1_year_location_summary.R`)
- 从 `PHENO.csv` 提取所有唯一的 year-location 组合
- 输出：`3_year_loc_summary/output/year_location_info.csv`（供后续模块使用）

**模块 4：天气数据下载** (`4_weather/code/1_get_met_apsim.R`)
- 为每个 year-location 下载 APSIM 格式的天气数据（`.met` 文件）
- 输出：`4_weather/output/[year-location].met`

**模块 5：APSIM 作物模型模拟** (`5_APSIM/code/1_apsim_sim.R`)
- 使用 APSIM 模型为每个 year-location 运行作物生长模拟
- 输出：`5_APSIM/output/simulations/[year-location].csv`（包含每日的生长阶段、温度、水分等变量）

**模块 6：环境协变量生成** (`6_ecov/code/1_get_env_cov.R`)
- 从 APSIM 模拟输出中提取不同生长阶段（如出苗期、开花期、灌浆期等）的环境变量
- 计算各阶段内的温度、降水、辐射等统计量（均值、总和、极值等）
- 输出：`6_ecov/output/ECOV.csv`（每行对应一个 year-location，列为各种环境协变量）

**模块 7：GDD 性状计算** (`7_GDD/code/1_add_GDD_to_pheno.R`)
- 基于 APSIM 输出的热时间（Thermal Time, TT）数据，计算开花相关性状的 GDD（Growing Degree Days）
- 为 `anthesis`（抽雄期）和 `silking`（吐丝期）计算 `anthesis_GDD` 和 `silking_GDD`
- 更新 `PHENO.csv`，添加 GDD 性状列

#### 运行 Pipeline

1. **准备原始数据**：从 G2F 网站下载数据，按上述结构放入 `source/` 目录

2. **按顺序运行各模块**（在 R 环境中，或使用 SLURM 提交作业）：
   ```r
   # 设置工作目录
   setwd("data_curation_and_ecov")
   
   # 模块 1：表型数据
   source("1_phenotypes/code/1_get_phenotype.R")
   
   # 模块 2：基因型数据（需要先运行 1_SNP_filter_codify.sub，再运行 2_LD_prune_genotypes.R）
   # 注意：1_SNP_filter_codify.sub 可能需要在高性能计算集群上提交
   source("2_genotypes/code/2_LD_prune_genotypes.R")
   
   # 模块 3：年份-地点摘要
   source("3_year_loc_summary/code/1_year_location_summary.R")
   
   # 模块 4：天气数据
   source("4_weather/code/1_get_met_apsim.R")
   
   # 模块 5：APSIM 模拟（可能需要较长时间）
   source("5_APSIM/code/1_apsim_sim.R")
   
   # 模块 6：环境协变量
   source("6_ecov/code/1_get_env_cov.R")
   
   # 模块 7：GDD 性状
   source("7_GDD/code/1_add_GDD_to_pheno.R")
   ```

3. **最终输出文件**：
   - `1_phenotypes/output/PHENO.csv`（包含 GDD 性状）
   - `2_genotypes/output/GENO.csv`
   - `2_genotypes/output/MAP.csv`
   - `6_ecov/output/ECOV.csv`

> **注意**：
> - 各模块脚本中的路径（如 `setwd(...)`）需要根据你的实际环境修改

---

### 1. 数据对齐

运行 `data_process/data_aligned.py` 脚本可以将**原始表型（PHENO）、基因型（GENO）、环境（ECOV）数据按样本一一对齐**，并生成后续训练使用的 `aligned_*.csv / .parquet` 文件。该脚本的核心流程如下：

- **基础路径与原始文件**
  - 在脚本开头通过 `os.chdir(...)` 和 `data_folder / output_folder` 指定原始数据与输出目录（请根据自己的数据路径进行修改）；
  - 读取表型数据与基因型数据 

- **GENO 与 PHENO 对齐（按基因型 ID）**

- **ECOV 与 PHENO 对齐（按 year-location）**

- **结果保存**

### 2. 输入数据文件约定

`main_gennet.py` 默认从数据目录读取以下文件（CSV 或 Parquet 二选一，若存在 Parquet 优先加载），使用时按需更改为实际目录：

- **基因型（SNP）数据示例**
  - CSV：`data/aligned_GENO.csv`
  - Parquet：`data/aligned_GENO.parquet`
  - 行 = 样本，列 = SNP 位点（第一列通常为样本 ID，将在预处理时被丢弃）

- **环境（ECOV）数据示例**
  - CSV：`data/aligned_ECOV.csv`
  - Parquet：`data/aligned_ECOV.parquet`
  - 行 = 样本，列 = 环境特征（第一列通常为样本 ID，将在预处理时被丢弃）

- **表型数据示例**
  - CSV：`data/South_PHENO.csv`
  - Parquet：`data/South_PHENO.parquet`
  - 至少包含以下目标性状列（顺序固定，对应多任务输出顺序）：
    - `yield`
    - `silking`
    - `anthesis`
    - `plants_stand`
    - `grain_moisture`
    - `anthesis_GDD`
    - `silking_GDD`

请确保：  
- 三个数据表在行上**已对齐同一批样本**（可通过运行 `data_process/data_aligned.py` 生成 `aligned_*.csv`，或自行实现等价逻辑）；  
- 行数一致，且上述 7 个表型列在表型数据中都存在。

---

## 完整运行示例流程

以下是从零开始到成功运行模型的完整步骤。

### 步骤 1：安装 R 语言及必需包

#### 1.1 安装 R（如果尚未安装）

**Windows:**
- 下载并安装 R 4.4.2

**Linux:**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install r-base r-base-dev

# 或使用 conda
conda install -c conda-forge r-base=4.4.2
```

#### 1.2 安装 R 必需包

打开 R 控制台或 RStudio，运行以下命令安装所需包：

```r
# 安装基础数据处理包
install.packages(c(
  "tidyverse",      # 数据处理
  "data.table",      # 高效数据读取
  "lubridate",      # 日期时间处理
  "jsonlite"        # JSON 处理
))

# 安装 APSIM 相关包（如果需要进行 APSIM 模拟）
install.packages(c(
  "apsimx",         # APSIM 模型接口
  "soilDB",         # 土壤数据库
  "spData",         # 空间数据
  "sf"              # 空间数据处理
))

# 安装可视化包（可选）
install.packages("ggrepel")
```

### 步骤 2：准备原始数据（如果从 G2F 数据开始）

#### 2.1 下载 G2F 原始数据

从 G2F 网站（https://www.genomes2fields.org/resources/）下载以下数据，并按结构放置：

#### 2.2 运行 R 数据处理 Pipeline

在 R 中按顺序运行各模块（或使用命令行）：

```r
# 设置工作目录
setwd("data_curation_and_ecov")

# 模块 1：表型数据整理
source("1_phenotypes/code/1_get_phenotype.R")
# 输出：1_phenotypes/output/PHENO.csv

# 模块 2：基因型数据处理
source("2_genotypes/code/2_LD_prune_genotypes.R")
# 输出：2_genotypes/output/GENO.csv, MAP.csv

# 模块 3：年份-地点摘要
source("3_year_loc_summary/code/1_year_location_summary.R")

# 模块 4：天气数据（如果需要）
source("4_weather/code/1_get_met_apsim.R")

# 模块 5：APSIM 模拟（如果需要，耗时较长）
source("5_APSIM/code/1_apsim_sim.R")

# 模块 6：环境协变量
source("6_ecov/code/1_get_env_cov.R")
# 输出：6_ecov/output/ECOV.csv

# 模块 7：GDD 性状计算
source("7_GDD/code/1_add_GDD_to_pheno.R")
# 更新：1_phenotypes/output/PHENO.csv（添加 GDD 列）
```

> **提示**：如果已有处理好的 `PHENO.csv`、`GENO.csv`、`ECOV.csv`，可跳过步骤 2，直接进入步骤 3。

### 步骤 3：数据对齐（Python）

#### 3.1 准备 Python 环境

```bash
# 创建虚拟环境
conda create -n GxE_env python==3.12.0
conda activate GxE_env
```

使用 requirements.txt 安装所有依赖：

```bash
pip install -r requirements.txt
```

#### 3.2 运行数据对齐脚本

修改 `data_process/data_aligned.py` 中的路径设置，然后运行：

```bash
python data_process/data_aligned.py
```

这将生成对齐后的文件：
- `aligned_GENO.csv` / `aligned_GENO.parquet`
- `aligned_ECOV.csv` / `aligned_ECOV.parquet`

#### 3.3 将文件移动到 data 目录

```bash
# 将生成的对齐文件复制到 data/ 目录
cp aligned_GENO.csv data/
cp aligned_ECOV.csv data/
cp [your_pheno_file].csv data/South_PHENO.csv  # 根据实际情况调整文件名
```

### 步骤 4：运行模型训练

#### 4.1 基础训练（不使用 GenNet 拓扑）

```bash
# 确保在项目根目录，且虚拟环境已激活
python main_gennet.py --epochs 100
```

这将：
- 自动加载 `data/` 目录下的数据
- 训练 100 个 epoch
- 在 `results/` 目录保存训练曲线和任务关系图
- 保存模型为 `final_model_gennet_YYYYMMDD_HHMMSS.pth`

### 步骤 5：查看结果

训练完成后，检查以下输出：

```bash
# 查看训练结果图像
ls results/
# 应包含：
# - train_loss_gennet_*.png
# - trait_*_metrics_gennet_*.png
# - task_relationships_gennet_*.png
# - task_weights_gennet_*.png

# 查看保存的模型
ls final_model_gennet_*.pth
```

</div>

<div id="lang-en" style="display: none;">

## Project Overview: Pheno Fusion GenNet Multitask Model

PyTorch-based multitask interaction model with optional GenNet topology (SNP→Gene) for **genotype × environment × multiple phenotypes**. Canonical entry: `main_gennet.py`.

- Multitask traits: `yield, silking, anthesis, plants_stand, grain_moisture, anthesis_GDD, silking_GDD`
- Optional GenNet topology via `--use_gennet` + `--topology_file`
- Task-relationship learning with visualized weights/relationships

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

</div>
