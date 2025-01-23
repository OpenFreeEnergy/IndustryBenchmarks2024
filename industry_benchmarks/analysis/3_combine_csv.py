import pandas as pd
import pathlib

csv_names = [
    "pymbar3_edge_data.csv",
    "pymbar3_calculated_dg_data.csv",
    "pymbar3_cumulative_data.csv",
    "pymbar4_edge_data.csv"
]

# map the folder to the dataset
folder_to_dataset = {}

for csv_type in csv_names:
    files = []
    for f in pathlib.Path("public_datasets_processed").glob(f"**/{csv_type}"):
        temp_df = pd.read_csv(f, index_col=0)
        # for csv with missing dataset name map the folder to the name of the dataset
        if csv_type == "pymbar3_edge_data.csv":
            folder_to_dataset[f.parent] = (temp_df.iloc[-1]["partner_id"], temp_df.iloc[-1]["dataset_name"])

        p_id, dataset_name = folder_to_dataset[f.parent]
        if "partner_id" not in temp_df.columns:
            temp_df["partner_id"] = p_id
        if "dataset_name" not in temp_df.columns:
            temp_df["dataset_name"] = dataset_name
        files.append(temp_df)

    master_file = pd.concat(files, ignore_index=True)
    master_file.sort_values(inplace=True, by="partner_id")
    master_file.reset_index(inplace=True, drop=True)
    master_file.to_csv(f"combined_{csv_type}", index=False)