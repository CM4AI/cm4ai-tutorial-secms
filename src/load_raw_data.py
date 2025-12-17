import os
import glob
import re
from pathlib import Path
import pandas as pd


def clean_input_labels(df_raw):
    cols = df_raw.columns.astype(str).tolist()
    cols = [re.sub(r"^PG\.", "", c) for c in cols]
    biosep_cols = [c for c in cols if "biosep" in c.lower()]
    rename_map = {c: f"biosep_{i+1}_quantity" for i, c in enumerate(biosep_cols)}
    cols = [rename_map.get(c, c) for c in cols]
    df_raw.columns = cols
    return df_raw


def get_disambiguated_dfs(df):
    # Remove rows with ambiguous gene/protein identification
    ambig_mask = (
            df["Genes"].str.contains(";", na=False) |
            df["ProteinAccessions"].str.contains(";", na=False) |
            df["ProteinNames"].str.contains(";", na=False)  
        )

    secms_data_clean = df.loc[~ambig_mask]
    secms_data_ambig = df.loc[ambig_mask]

    return secms_data_clean, secms_data_ambig


def load_secms_data(base_path):
    filename_pattern = r"[_-](?P<category>[A-Za-z]+)_(?P<expnum>\d+)_Report"

    data_sets = [p.name for p in Path(base_path).iterdir() if p.is_dir()]
    
    dfs = []
    for data_set in data_sets:
        data_path = os.path.join(base_path, data_set)
        cell_line = data_set.split("_")[-1]

        data_files = glob.glob(os.path.join(data_path, '*.tsv'))

        for data_file in data_files:
            filename = Path(data_file).stem
            m = re.search(filename_pattern, filename)

            category = None
            expnum = None
            if m:
                category = m.group("category")
                expnum = int(m.group("expnum"))

            df = pd.read_csv(data_file, sep="\t")
            df = clean_input_labels(df)
            df["category"] = category
            df["expnum"] = expnum
            df["cell_line"] = cell_line
            df["source_file"] = data_file

            dfs.append(df)

    secms_data = pd.concat(dfs, ignore_index=True, sort=False)

    return secms_data