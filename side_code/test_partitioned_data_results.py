import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from scipy.stats import chisquare


def evaluate_data_per_gene(data):
    for index, row in data.iterrows():
        if row.get("full_data_counts") and row["actual_sample_pct"] == 0.05:
            dataset = re.findall('\_([0-9A-Za-z]+)\.((fasta)|(phy))', row["curr_msa_version_folder"])[0][0]
            if dataset in ["NagyA1", "YangA8", "PrumD6"]:
                    print(dataset)
                #try:
                    print(row["lasso_test_R^2"])
                    print(row["actual_sample_pct"])
                    full = [int(n) for n in row["full_data_counts"][1:-1].split(" ") if len(n) > 0][1:]
                    # print(full)
                    chosen_expected = [int(n) * row["actual_sample_pct"] for n in full]
                    # print(chosen_expected)
                    not_chosen_expected = [int(n) * (1 - row["actual_sample_pct"]) for n in full]
                    # print(not_chosen_expected)
                    mean_rate4site = [float(val) for val in row["partition_mean_rates"][1:-1].split(',')[1:]]
                    print("n_categories=",len(mean_rate4site))
                    # mean_coeff = [float(val) for val in row["partition_mean_coeff"][1:-1].split(',')[1:]]
                    chosen = [int(n) for n in row["observed_partition_counts"][1:-1].split(" ") if len(n) > 0][1:]
                    # print(chosen)
                    not_chosen = [f - n for f, n in zip(full, chosen)]
                    # print(not_chosen)
                    chi_square_statistics = chisquare(f_obs=chosen + not_chosen,
                                                      f_exp=chosen_expected + not_chosen_expected)
                    print(chi_square_statistics)
                    # df = pd.DataFrame({"chosen_expected": chosen_expected, "observed": chosen},index = range(len(full)))
                    # df.plot.bar(rot=0)
                    # plt.show()
                    chosen_by_expected = np.divide(np.array(chosen), np.array(chosen_expected))
                    print(np.corrcoef(x=np.array(chosen_by_expected), y=np.array(mean_rate4site)))
                    df2 = pd.DataFrame({"chosen_by_expected": chosen_by_expected, "rate4site": mean_rate4site},
                                       index=range(len(full)))
                    # df2.plot.scatter(x="rate4site", y = "chosen_by_expected", rot=0)
                    # plt.show()

                #except:
                #    print(f"Could not be applied to dataset {dataset}")


def main():
    print("Not partitioned_data:")
    no_partition_data_aa = pd.read_csv(
        "/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_no_partition.tsv", sep='\t')
    no_partition_data_dna = pd.read_csv(
        "/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_no_partition_dna.tsv", sep='\t')
    no_partition_combined = pd.concat([no_partition_data_aa, no_partition_data_dna], ignore_index=True)
    evaluate_data_per_gene(no_partition_combined)
    print("\n\npartitioned_data:")
    partition_data_aa = pd.read_csv(
        "/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_partition.tsv", sep='\t')
    partition_data_dna = pd.read_csv(
        "/Users/noa/Workspace/lasso_positions_sampling_results/revision_files/test_partition_dna.tsv", sep='\t')
    partition_combined = pd.concat([partition_data_aa, partition_data_dna], ignore_index=True)
    evaluate_data_per_gene(partition_combined )


if __name__ == "__main__":
    main()
