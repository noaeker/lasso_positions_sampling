
import pickle
from config import SEED
from help_functions import create_or_clean_dir, create_dir_if_not_exists, update_dict_with_a_suffix
import os
import shutil
import logging
from raxml import generate_n_random_tree_topology_constant_brlen,RF_between_two_newick
from generate_SPR import SPR_analysis
from config import IGNORE_COLS_IN_CSV

def generate_or_copy_random_starting_tree(i, curr_run_directory, curr_msa_stats):
    seed = SEED + i+ curr_msa_stats["start_of_starting_tree_ind"]
    random_tree_folder = os.path.join(curr_run_directory, "RANDOM_starting_tree_" + str(i))
    create_or_clean_dir(random_tree_folder)
    starting_tree_path = os.path.join(random_tree_folder, ".raxml.startTree")
    baseline_starting_tree_path = starting_tree_path.replace(curr_msa_stats["run_prefix"],
                                                             curr_msa_stats["spr_baseline_run_prefix"])
    if os.path.exists(baseline_starting_tree_path):
        shutil.copyfile(baseline_starting_tree_path, starting_tree_path)
    if not os.path.exists(starting_tree_path):
        starting_tree_path, elapsed_running_time = generate_n_random_tree_topology_constant_brlen(
            curr_msa_stats=curr_msa_stats, n=1, alpha=curr_msa_stats["alpha"],
            original_file_path=curr_msa_stats["local_alignment_path"],
            curr_run_directory=random_tree_folder, seed=seed)
    return starting_tree_path



def lasso_spr_pipeline_per_training_data(curr_msa_stats, training_size, brlen_run_directory, starting_tree_path, lasso_configurations_per_training_size, brlen_generator_name):
    curr_msa_stats["actucal_training_size"] = training_size
    curr_training_size_directory = os.path.join(brlen_run_directory, str(training_size))
    create_dir_if_not_exists(curr_training_size_directory)
    lasso_config_per_brlen_and_t_size = lasso_configurations_per_training_size[brlen_generator_name][
        training_size]
    curr_msa_stats["curr_training_size_lasso_configurations"] = lasso_config_per_brlen_and_t_size
    logging.info(f"    ****Using training size {training_size}\n")
    spr_search_configurations = []

    lasso_thresholds_during_search = [float(t) for t in curr_msa_stats['lasso_thresholds_search'].split("_")]
    top_ind_to_test_per_phase = [int(t) for t in curr_msa_stats['top_ind_to_test_per_phase'].split("_")]

    for i, per_phase_search_data in enumerate(zip(lasso_thresholds_during_search,top_ind_to_test_per_phase)):
        threshold, top_ind = per_phase_search_data
        lasso_results = lasso_config_per_brlen_and_t_size[float(threshold)]
        lasso_results["top_ind"] = top_ind
        update_dict_with_a_suffix(curr_msa_stats, lasso_results, suffix=f"_phase_{i}")
        spr_search_configurations.append(lasso_results)
    lasso_based_spr_results = SPR_analysis(
        spr_search_configurations,starting_tree_path, curr_msa_stats,
        curr_run_directory=curr_training_size_directory, full_run=False)

    lasso_based_spr_results_print = {k: lasso_based_spr_results[k] for k in lasso_based_spr_results.keys() if all(x not in k for x in ["path", "newick","pval"])}




    logging.info(f"     *****Lasso results were succesfully obtained: \n{lasso_based_spr_results_print}\n ")
    curr_msa_stats.update(lasso_based_spr_results)
    rf_folder = os.path.join(curr_training_size_directory, "rf_calculations")
    create_dir_if_not_exists(rf_folder)
    curr_msa_stats["rf_best_naive_vs_best_lasso"] = RF_between_two_newick(rf_folder, "best_vs_lasso",curr_msa_stats["naive_SPR_best_tree_newick"], curr_msa_stats["final_phase_best_tree_newick"])


def perform_spr_pipeline(training_size_options, brlen_generators, curr_msa_stats,
                         lasso_configurations_per_training_size,
                         job_csv_path, all_msa_results):
    spr_searches_run_directory = os.path.join(curr_msa_stats["curr_msa_version_folder"], "spr_results")
    create_dir_if_not_exists(spr_searches_run_directory)
    logging.info("\n\nStarting SPR searches")
    n_starting_trees = curr_msa_stats["n_random_starting_trees"] if not curr_msa_stats[
        "use_spr_parsimony_starting_tree"] else 1
    for i in range(n_starting_trees):
        logging.info(f" *Starting tree {i}:")
        starting_tree_run_directory = os.path.join(spr_searches_run_directory,
                                                   'starting_tree_{i}'.format(i=i))
        create_dir_if_not_exists(starting_tree_run_directory)
        if curr_msa_stats["use_spr_parsimony_starting_tree"]:
            logging.info(f"Using a parsimony starting tree to spr search")
            starting_tree_path = curr_msa_stats["parsimony_optimized_tree_path"]
        else:
            starting_tree_path = generate_or_copy_random_starting_tree(i, starting_tree_run_directory, curr_msa_stats)
        curr_msa_stats["starting_tree_path"] = starting_tree_path
        naive_spr_results_dump = os.path.join(starting_tree_run_directory, 'naive_spr.dump')
        naive_spr_results_dump_baseline = naive_spr_results_dump.replace(curr_msa_stats["run_prefix"],
                                                                         curr_msa_stats["spr_baseline_run_prefix"])
        if os.path.exists(naive_spr_results_dump_baseline):
            logging.info("**Using full data naive SPR dump in {}".format(naive_spr_results_dump_baseline))
            with open(naive_spr_results_dump_baseline, 'rb') as handle:
                naive_spr_results = pickle.load(handle)
                logging.info(f"Naive SPR extracted from dump are:\n{naive_spr_results}\n")
        else:
            logging.info("  **Running full data naive SPR from beggining")
            curr_starting_tree_full_run_directory = os.path.join(starting_tree_run_directory, "spr_full_data_results")
            create_dir_if_not_exists(curr_starting_tree_full_run_directory)
            naive_spr_results = SPR_analysis(None,starting_tree_path,
                                             curr_msa_stats,
                                             curr_run_directory=curr_starting_tree_full_run_directory,
                                             full_run=True)
            with open(naive_spr_results_dump, 'wb') as handle:
                pickle.dump(naive_spr_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        curr_msa_stats.update(naive_spr_results)
        logging.info("  **Running Lasso-based SPR searches using the thresholds:")
        for brlen_generator_name in brlen_generators:
            brlen_run_directory = os.path.join(starting_tree_run_directory, brlen_generator_name)
            create_dir_if_not_exists(brlen_run_directory)
            curr_msa_stats["brlen_generator"] = brlen_generator_name
            logging.info(f"   ***Using brlen_generator {brlen_generator_name} ")
            for training_size in training_size_options:

                lasso_spr_pipeline_per_training_data(curr_msa_stats, training_size, brlen_run_directory,
                                                      starting_tree_path,
                                                     lasso_configurations_per_training_size, brlen_generator_name)
                all_msa_results = all_msa_results.append({k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                                          k not in IGNORE_COLS_IN_CSV
                                                          }, ignore_index=True)
                all_msa_results.to_csv(job_csv_path,index = False,sep ='\t')
    return all_msa_results
