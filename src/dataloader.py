import yaml
import pandas as pd

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def load_dataset(report_path):
    df = pd.read_csv(
    report_path,
    sep="\t",
    low_memory=False,
    dtype=str,
    encoding="utf-8",      # 如有乱码可尝试 "latin1"
    )
    return df