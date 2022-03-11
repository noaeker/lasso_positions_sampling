import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.stats import chisquare

no_partition_data = pd.read_csv("/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_no_partition.tsv", sep = '\t')
#partition_data = pd.read_csv("/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_partition.tsv",sep = '\t')

for index,row in no_partition_data.iterrows():
    try:
        dataset = re.findall('\_([0-9A-Za-z]+)\.((fasta)|(phy))',row["curr_msa_version_folder"])[0][0]
        print(dataset)
        print(row["actual_sample_pct"])
        full = [int(n) for n in row["full_data_counts"][1:-1].split(" ") if len(n)>0][1:]
        chosen_expected = [int(n)*row["actual_sample_pct"]  for n in full]
        not_chosen_expected = [int(n)*(1-row["actual_sample_pct"])  for n in full]
        mean_rate4site = [float(val) for val in row["partition_mean_rates"][1:-1].split(',')[1:]]
        #mean_coeff = [float(val) for val in row["partition_mean_coeff"][1:-1].split(',')[1:]]
        chosen = [int(n) for n in row["observed_partition_counts"][1:-1].split(" ") if len(n) > 0][1:]
        not_chosen = [f-n for f,n in zip(full,chosen)]
        chi_square_statistics = chisquare(f_obs= chosen+not_chosen, f_exp=chosen_expected+not_chosen_expected)
        #df = pd.DataFrame({"chosen_expected": chosen_expected, "observed": chosen},index = range(len(full)))
        #df.plot.bar(rot=0)
        #plt.show()
        chosen_by_expected = np.divide(np.array(chosen), np.array(chosen_expected))
        print(np.corrcoef(x=np.array(chosen_by_expected), y=np.array(mean_rate4site)))
        df2 = pd.DataFrame({"chosen_by_expected": chosen_by_expected, "rate4site": mean_rate4site}, index=range(len(full)))
        df2.plot.scatter(x="rate4site", y = "chosen_by_expected", rot=0)
        plt.show()
        #df2.plot.scatter(x="coeff", y="chosen_by_expected", rot=0)
        #plt.show()
    except:
        print(f"Could not be applied to dataset {dataset}")
