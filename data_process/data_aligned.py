import pandas as pd
import os

os.chdir("/mnt/uict_data/wanghr/Maize_Hub_py")

data_folder = "data_ori"
output_folder = "data_processed"

def load_or_convert_data(csv_file, parquet_file):
    """Load data from Parquet if exists, otherwise read from CSV and convert to Parquet."""
    try:
        if os.path.exists(parquet_file):
            print(f"Loading GENO from {parquet_file}...")
            data = pd.read_parquet(parquet_file)
            print("Loaded GENO from parquet success.")
        else:
            print(f"Parquet file not found, reading GENO from {csv_file} and converting to parquet...")
            data = pd.read_csv(csv_file)
            print("Read GENO from CSV success.")

            # 将数据保存为Parquet格式
            data.to_parquet(parquet_file)
            print(f"Converted and saved GENO to {parquet_file}.")

        # 检查数据是否为空
        if data.empty:
            print("Warning: Loaded data is empty.")
            return None

        return data

    except Exception as e:
        print(f"Error loading or converting data: {e}")
        raise

# 读取数据
pheno_df = pd.read_csv("data_processed/North_PHENO.csv")
print("read South_PHENO success")
csv_file_path = "data_processed/aligned_North_GENO.csv"
parquet_file_path = "data_processed/aligned_North_GENO.parquet"
geno_df = load_or_convert_data(csv_file_path, parquet_file_path)
dfcolumns = geno_df.columns
print("read GENO success")
ecov_df = pd.read_csv("data_ori/ECOV.csv")
print("read ECOV success")

# 设置GENO的索引为genotype列（假设GENO中存在该列）
geno_df.set_index(geno_df.columns[0], inplace=True)

# 检查GENO中是否存在PHENO的所有genotype
missing_geno = pheno_df[~pheno_df["genotype"].isin(geno_df.index)]
if not missing_geno.empty:
    print("警告：以下genotype在GENO中不存在：\n", missing_geno)

# 提取与PHENO对应的行（保持PHENO的顺序）
aligned_geno = geno_df.loc[pheno_df["genotype"]]
aligned_geno.to_csv("aligned_GENO.csv", index=True)  # 保留genotype作为索引
print("aligned_geno.to_csv success")

# 2. 处理ECOV文件（根据PHENO的year和location列匹配）
# 在PHENO中创建year-location列
pheno_df["year_location"] = pheno_df["year"].astype(str) + "-" + pheno_df["location"].astype(str)
ecov_df.set_index(ecov_df.columns[0], inplace=True)

# 提取与PHENO对应的行（保持PHENO的顺序）
aligned_ecov = ecov_df.loc[pheno_df["year_location"]]
aligned_ecov.to_csv("aligned_North_ECOV.csv", index=True)  # 保留year_location作为索引
print("aligned_ecov.to_csv success")

print("对齐完成！已生成aligned_GENO.csv和aligned_ECOV.csv")