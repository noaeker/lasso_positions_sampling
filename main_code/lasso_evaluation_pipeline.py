
from training_and_test_set_generation import Lasso_training_and_test
import logging
from config import IGNORE_COLS_IN_CSV
from lasso_model_analysis import compare_lasso_to_naive_approaches_on_test_set
from partitioned_analysis import get_mean_param_per_group
import matplotlib.pyplot  as plt
import pickle
from scipy.stats import chisquare
import os
import numpy as np




def get_lasso_configurations(curr_msa_version_folder, args, brlen_generators, curr_msa_stats, training_size_options):
    curr_msa_version_lasso_dump = os.path.join(curr_msa_version_folder, 'lasso.dump')
    curr_msa_version_lasso_dump_baseline = curr_msa_version_lasso_dump.replace(args.run_prefix,
                                                                               args.lasso_baseline_run_prefix)
    if os.path.exists(curr_msa_version_lasso_dump_baseline):
        with open(curr_msa_version_lasso_dump_baseline, 'rb') as handle:
            logging.info(
                "Using lasso dump files in {} ".format(curr_msa_version_lasso_dump))
            lasso_configurations_per_training_size = pickle.load(handle)
    else:
        lasso_configurations_per_training_size = Lasso_training_and_test(brlen_generators, curr_msa_stats,
                                                                         training_size_options,
                                                                         args.random_trees_test_size)
        with open(curr_msa_version_lasso_dump, 'wb') as handle:
            pickle.dump(lasso_configurations_per_training_size, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return lasso_configurations_per_training_size

def perform_only_lasso_pipeline(training_size_options, brlen_generators, curr_msa_stats,
                                lasso_configurations_per_training_size,
                                job_csv_path, all_msa_results, curr_run_directory):
    for brlen_generator_name in brlen_generators:
        curr_msa_stats["brlen_generator"] = brlen_generator_name
        for training_size in training_size_options:
            curr_msa_stats["actual_training_size"] = training_size
            lasso_results = lasso_configurations_per_training_size[brlen_generator_name][training_size]
            for threshold in lasso_results:
                curr_msa_stats["actual_sample_pct"] = threshold
                curr_msa_stats.update(lasso_results[threshold])
                logging.info(
                    "only evaluating lasso on brlen {} and training size {}: ".format(brlen_generator_name,
                                                                                      training_size))
                lasso_evaluation_result = {k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                           k not in IGNORE_COLS_IN_CSV
                                           }
                if curr_msa_stats["compare_lasso_to_naive"]:
                    lasso_comparisons_results = compare_lasso_to_naive_approaches_on_test_set(curr_msa_stats, curr_run_directory, threshold)
                    lasso_evaluation_result.update(lasso_comparisons_results)
                if curr_msa_stats["compare_loci_gene_distribution"] and curr_msa_stats['per_loci_partition'] is not None:
                    chosen_locis = curr_msa_stats["lasso_chosen_locis"]
                    logging.info(f"Partition count = {curr_msa_stats['per_loci_partition']}")
                    partitions_count_arr = np.bincount(curr_msa_stats["per_loci_partition"])
                    partitioned_rates = get_mean_param_per_group(curr_msa_stats["rate4site_scores"],curr_msa_stats['per_loci_partition'])
                    partitioned_coeffs = get_mean_param_per_group(curr_msa_stats["lasso_chosen_weights"],
                                                                 np.take(np.array(curr_msa_stats['per_loci_partition']),curr_msa_stats["lasso_chosen_locis"]))
                    expected_chosen_locis_count_arr = (np.bincount(curr_msa_stats["per_loci_partition"])*threshold).astype(int)
                    chosen_locis_partitions_count_arr = np.bincount(curr_msa_stats["per_loci_partition"][chosen_locis])
                    chosen_locis_partitions_count_arr = np.pad(chosen_locis_partitions_count_arr,(0,len(partitions_count_arr)-len(chosen_locis_partitions_count_arr)), mode = 'constant')
                    obs_vs_expected = np.divide(chosen_locis_partitions_count_arr,partitions_count_arr)
                    rates = [partitioned_rates.get(i)  for i in range(len(obs_vs_expected))]
                    try:
                        chi_square_statistics = chisquare(chosen_locis_partitions_count_arr,f_exp = expected_chosen_locis_count_arr )

                    except:
                        chi_square_statistics = -1
                    lasso_evaluation_result["expected_partition_counts"] = expected_chosen_locis_count_arr
                    lasso_evaluation_result["partition_mean_rates"] = rates
                    lasso_evaluation_result["partition_mean_coeff"] =  partitioned_coeffs
                    lasso_evaluation_result["full_data_counts"] = partitions_count_arr
                    lasso_evaluation_result["observed_partition_counts"] = chosen_locis_partitions_count_arr
                    lasso_evaluation_result["chi_square_partition"] = chi_square_statistics


                all_msa_results = all_msa_results.append(lasso_evaluation_result, ignore_index=True)
                all_msa_results.to_csv(job_csv_path, index=False,sep ='\t')
    return all_msa_results
