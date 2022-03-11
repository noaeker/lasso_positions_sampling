import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.stats import chisquare

no_partition_data = pd.read_csv("/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_no_partition.tsv", sep = '\t')
partition_data = pd.read_csv("/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_partition.tsv",sep = '\t')

for index,row in no_partition_data.iterrows():
        dataset = re.findall('\_([0-9A-Za-z]+)\.((fasta)|(phy))',row["curr_msa_version_folder"])
        try:
            print(dataset)
            print(row["actual_sample_pct"])
            full = [int(n) for n in row["full_data_counts"][1:-1].split(" ") if len(n)>0][1:]
            chosen_expected = [int(n)*row["actual_sample_pct"]  for n in full]
            not_chosen_expected = [int(n)*(1-row["actual_sample_pct"])  for n in full]
            chosen = [int(n) for n in row["observed_partition_counts"][1:-1].split(" ") if len(n) > 0][1:]
            not_chosen = [f-n for f,n in zip(full,chosen)]
            print(chosen)
            print(chosen_expected)
            chi_square_statistics = chisquare(f_obs= chosen+not_chosen, f_exp=chosen_expected+not_chosen_expected)
            print(chi_square_statistics)
            #df = pd.DataFrame({"expected": expected, "observed": observed},index = range(len(full)))
            #df .plot.bar(rot=0)
            #plt.show()
        except:
            f"Could not be applied to {dataset}"
