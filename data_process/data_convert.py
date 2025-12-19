import pandas as pd
import os

# os.chdir("/mnt/uict_data/wanghr/Maize_DL_Server_16")

# data_folder = "data"

def load_or_convert_data(csv_file, parquet_file):
    """Load data from Parquet if exists, otherwise read from CSV and convert to Parquet."""
    try:
        if os.path.exists(parquet_file):
            print(f"Loading file from {parquet_file}...")
            data = pd.read_parquet(parquet_file)
            print("Loaded file from parquet success.")
        else:
            print(f"Parquet file not found, reading file from {csv_file} and converting to parquet...")
            data = pd.read_csv(csv_file)
            print("Read file from CSV success.")

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




# ecov_df = pd.read_csv("data/aligned_ECOV.csv")
# geno_df = pd.read_csv("data/aligned_GENO.csv")
#
# parquet_file_path_ecov = "data/aligned_ECOV.parquet"
# parquet_file_path_geno = "data/aligned_GENO.parquet"
# ecov_df.to_parquet(parquet_file_path_ecov)
# geno_df.to_parquet(parquet_file_path_geno)
# print(f"Converted and saved GENO to {parquet_file_path_geno}.")
#
# print("read success")