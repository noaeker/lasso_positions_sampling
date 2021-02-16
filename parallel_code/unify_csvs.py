import pandas as pd
import os
import logging

def add_csvs_content(csvs_path_list, unified_csv_path):
    existing_df = [pd.read_csv(unified_csv_path)] if os.path.exists(unified_csv_path) else []
    existing_df_size = pd.read_csv(unified_csv_path).size if os.path.exists(unified_csv_path) else 0
    logging.info('Existing df size is: {}'.format(existing_df_size))
    non_empty_df = [pd.read_csv(f) for f in csvs_path_list if not pd.read_csv(f).empty]
    combined_df = pd.concat(non_empty_df+ existing_df,sort=False)
    combined_df_size = combined_df.size
    logging.info('Combined df size is: {}'.format(combined_df_size))
    combined_df.to_csv(unified_csv_path,index=False)
    return combined_df


unified_csv_path = "/Users/noa/Workspace/lasso_positions_sampling_results/curr_spr_test.csv"
for i in range(50):
    combined_df = add_csvs_content(["/Users/noa/Workspace/lasso_positions_sampling_results/tmp_csvs/{}.csv".format(i)],unified_csv_path)
