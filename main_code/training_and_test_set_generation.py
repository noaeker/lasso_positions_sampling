#from lasso_model_analysis import *
import pickle
import os
from help_functions import create_dir_if_not_exists
import logging
import random
from config import SEED
import pandas as pd
import shutil
from raxml import raxml_optimize_trees_for_given_msa,raxml_compute_tree_per_site_ll, generate_n_random_tree_topology_constant_brlen
from lasso_model_analysis import apply_lasso_on_sitelh_data_and_update_statistics
from spr_prune_and_regraft import add_internal_names,get_possible_spr_moves, generate_neighbour,generate_tree_object_from_newick,assign_brlen_to_tree_object


def Lasso_test_set(curr_msa_stats, random_trees_test_size, Lasso_folder, random_trees_folder, test_seed):
    if not curr_msa_stats["no_test_set"]:
        test_random_trees_path, test_random_tree_generation_time = generate_n_random_topologies_constant_brlen(
            random_trees_test_size, random_trees_folder,
            curr_msa_stats,
            "test", seed=test_seed)
        test_optimization_folder = os.path.join(Lasso_folder,
                                                "test_{}_random_trees_eval".format(random_trees_test_size))
        create_dir_if_not_exists(test_optimization_folder)
        test_ll_values, optimized_test_topologies_path = generate_optimized_tree_topologies_for_testing(curr_msa_stats,
                                                                                                        test_random_trees_path,
                                                                                                        test_optimization_folder)
        curr_msa_stats["test_ll_values"] = test_ll_values
        curr_msa_stats["optimized_test_topologies_path"] = optimized_test_topologies_path
    else:
        test_ll_values, optimized_test_topologies_path = None, None
    return test_ll_values, optimized_test_topologies_path


def generate_specific_brlen_training_set(brlen_generator_name, Lasso_folder, brlen_generators, training_size,
                                         curr_msa_stats, training_random_trees_path):
    brlen_run_directory = os.path.join(Lasso_folder, brlen_generator_name)
    create_dir_if_not_exists(brlen_run_directory)
    brlen_generator_func = brlen_generators.get(brlen_generator_name)
    training_size_directory = os.path.join(brlen_run_directory,
                                           "training_{}_random_tree_eval".format(training_size))
    create_dir_if_not_exists(training_size_directory)
    training_output_csv_path = os.path.join(training_size_directory,
                                            "training" + ".csv")
    training_dump = os.path.join(training_size_directory, 'training_set.dump')
    training_dump_baseline = training_dump.replace(curr_msa_stats["run_prefix"],
                                                   curr_msa_stats["training_set_baseline_run_prefix"])
    alternative_data_path = (
        os.path.join(curr_msa_stats["curr_msa_version_folder"], "actual_search_training_df.csv")).replace(
        curr_msa_stats["run_prefix"], curr_msa_stats["alternative_training_prefix"])

    if os.path.exists(training_dump_baseline):
        logging.info("Using trainng results in {}".format(training_dump_baseline))
        with open(training_dump_baseline, 'rb') as handle:
            training_results = pickle.load(handle)
            training_sitelh, training_eval_time = training_results["training_sitelh"], training_results[
                "training_eval_time"]
    elif (os.path.exists(alternative_data_path) and curr_msa_stats["alternative_training"]):
        logging.info(f"Using alternative training set in {alternative_data_path}")
        training_sitelh = pd.read_csv(alternative_data_path)
        training_eval_time = 0
        random.seed(SEED)
        training_sitelh = training_sitelh.sample(frac=1)
    else:
        training_sitelh, training_eval_time = generate_per_site_ll_on_random_trees_for_training(
            curr_msa_stats=curr_msa_stats,
            random_trees_path=training_random_trees_path,
            brlen_generator_func=brlen_generator_func,
            curr_run_directory=training_size_directory,
            output_csv_path=training_output_csv_path)
        training_results = {"training_sitelh": training_sitelh, "training_eval_time": training_eval_time}
        with open(training_dump, 'wb') as handle:
            pickle.dump(training_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logging.debug(
        "Done evaluating topologies based on {} branch lengths. It took {} seconds".format(brlen_generator_name,
                                                                                           training_eval_time))
    return training_sitelh, training_eval_time, training_size_directory


def Lasso_training_and_test(brlen_generators, curr_msa_stats, training_size_options, random_trees_test_size):
    Lasso_folder = os.path.join(curr_msa_stats["curr_msa_version_folder"], "Lasso_folder")
    create_dir_if_not_exists(Lasso_folder)
    logging.info("Generating Lasso folder in {}".format(Lasso_folder))
    random_trees_folder = os.path.join(Lasso_folder, "random_tree_generation")
    curr_msa_stats["Lasso_folder"] = Lasso_folder
    curr_msa_stats["random_trees_folder"] = random_trees_folder
    create_dir_if_not_exists(random_trees_folder)
    start_seed_random_trees = SEED
    test_ll_values, optimized_test_topologies_path = Lasso_test_set(curr_msa_stats, random_trees_test_size,
                                                                    Lasso_folder, random_trees_folder,
                                                                    start_seed_random_trees)
    max_training_size = max(training_size_options)
    if curr_msa_stats["use_spr_neighbours_training"]:
        training_random_trees_path, training_tree_generation_elapsed_running_time = generate_n_random_spr_neighbours_topologies_orig_brlen(
            max(training_size_options), random_trees_folder,
            curr_msa_stats, "training",
            seed=start_seed_random_trees)
    else:
        training_random_trees_path, training_tree_generation_elapsed_running_time = generate_n_random_topologies_constant_brlen(
            max(training_size_options), random_trees_folder,
            curr_msa_stats, "training",
            seed=start_seed_random_trees)

    run_configurations = {}
    logging.info(f'Generating Lasso results:')
    for brlen_generator_name in brlen_generators:
        logging.info(f' *Obtaining training data for brlen {brlen_generator_name}:')
        training_sitelh, training_eval_time, training_full_size_directory = generate_specific_brlen_training_set(
            brlen_generator_name, Lasso_folder, brlen_generators, max_training_size, curr_msa_stats,
            training_random_trees_path)

        for training_size in training_size_options:
            logging.info(
                f'  **Applying Lasso for various alphas on current trimmed training data of size: {training_size}')
            training_sitelh_trimmed = training_sitelh.iloc[:training_size].copy()
            trimmed_training_directory = os.path.join(training_full_size_directory, f"trimmed_{training_size}")
            create_dir_if_not_exists(trimmed_training_directory)
            Lasso_results = apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats,
                                                                             curr_run_directory=trimmed_training_directory,
                                                                             sitelh_training_df=training_sitelh_trimmed,
                                                                             test_optimized_trees_path=optimized_test_topologies_path)  # calculating positions_weight

            for threshold in Lasso_results:
                Lasso_results[threshold].update(
                    {'full_training_random_trees_generation_time': training_tree_generation_elapsed_running_time,
                     'full_size_training_evaluation_time': training_eval_time,
                     'lasso_training_size': training_size,
                     'lasso_brlen_generator': brlen_generator_name
                     })
            if brlen_generator_name not in run_configurations:
                run_configurations[brlen_generator_name] = {}
            run_configurations[brlen_generator_name][training_size] = Lasso_results
    shutil.rmtree(random_trees_folder)
    shutil.rmtree(training_full_size_directory)
    return run_configurations


def generate_n_random_spr_neighbours_topologies_orig_brlen(n, curr_run_directory, curr_msa_stats, name, seed):
    total_trees = 0
    all_tree_objects = []
    curr_seed = seed
    random_tree_path, elapsed_running_time = generate_n_random_topologies_constant_brlen(1, curr_run_directory,
                                                                                         curr_msa_stats, name,
                                                                                         seed=curr_seed)
    trees_ll_on_data, tree_object, elapsed_running_time = raxml_optimize_trees_for_given_msa(
        curr_msa_stats["local_alignment_path"], "opt_first_training_tree", random_tree_path, curr_msa_stats,
        curr_run_directory, opt_brlen=True, weights=None, return_trees_file=False,
        n_cpus=curr_msa_stats["n_cpus_training"])
    add_internal_names(tree_object)
    tree_object.get_tree_root().name = "ROOT"
    curr_tree_object= tree_object
    while total_trees<n:
        #print(f"total_trees = {total_trees}")
        spr_moves = get_possible_spr_moves(curr_tree_object,
                                                              rearr_dist=curr_msa_stats["rearr_dist"])
        random.seed(curr_seed)
        random_spr_move = random.choice(spr_moves)
        random_spr_neighbour = generate_neighbour(curr_tree_object, random_spr_move )
        all_tree_objects = all_tree_objects+ [random_spr_neighbour]
        total_trees = total_trees + 1
        curr_seed = curr_seed + 1
        add_internal_names(random_spr_neighbour)
        random_spr_neighbour.get_tree_root().name = "ROOT"
        curr_tree_object = random_spr_neighbour
    logging.info(f"Used {n} random optimized trees to generate training data!")
    trees_newick = "\n".join([tree.write(format=1) for tree in all_tree_objects[:n]])
    trees_eval_path = os.path.join(curr_run_directory, "all_training_neighbours_trees")
    with open(trees_eval_path, 'w') as EVAL_TREES:
        EVAL_TREES.write(trees_newick)
    return trees_eval_path, elapsed_running_time

def generate_n_random_topologies_constant_brlen(n, curr_run_directory, curr_msa_stats, name, seed):
    local_file_path = curr_msa_stats.get("local_alignment_path")
    basic_directory = os.path.join(curr_run_directory, "{name}_{n}".format(name=name, n=n))
    create_dir_if_not_exists(basic_directory)
    seed = seed
    alpha = curr_msa_stats["alpha"]
    random_tree_path, elapsed_running_time = generate_n_random_tree_topology_constant_brlen(n, alpha, local_file_path,
                                                                                             basic_directory,
                                                                                             curr_msa_stats, seed=seed)
    return random_tree_path, elapsed_running_time


def generate_optimized_tree_topologies_for_testing(curr_msa_stats, test_random_trees_path, curr_run_directory):
    local_file_path = curr_msa_stats.get("local_alignment_path")
    test_dump = os.path.join(curr_run_directory, 'test_set.dump')
    test_dump_baseline = test_dump.replace(curr_msa_stats["run_prefix"],
                                           curr_msa_stats["test_set_baseline_run_prefix"])

    if os.path.exists(test_dump_baseline):
        logging.info("Using test results in {}".format(test_dump_baseline))
        with open(test_dump_baseline, 'rb') as handle:
            test_results = pickle.load(handle)
            trees_ll_on_data, optimized_trees_path = test_results["test_ll_values"], test_results[
                "optimized_test_topologies_path"]
    else:
        trees_ll_on_data, optimized_trees_path, elapsed_time = raxml_optimize_trees_for_given_msa(local_file_path,
                                                                                                  "test_opt",
                                                                                                  test_random_trees_path,
                                                                                                  curr_msa_stats,
                                                                                                  curr_run_directory,
                                                                                                  weights=False,
                                                                                                  return_trees_file=True,
                                                                                                  n_cpus=curr_msa_stats[
                                                                                                      "n_cpus_training"])
        test_results = {"test_ll_values": trees_ll_on_data, "optimized_test_topologies_path": optimized_trees_path}
        with open(test_dump, 'wb') as handle:
            pickle.dump(test_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return trees_ll_on_data, optimized_trees_path


def generate_per_site_ll_on_random_trees_for_training(curr_msa_stats, random_trees_path, brlen_generator_func,
                                                      curr_run_directory, output_csv_path):
    raxml_ll_eval_directory = os.path.join(curr_run_directory, "raxml_training_per_site_ll_eval")
    create_dir_if_not_exists(raxml_ll_eval_directory)
    local_file_path = curr_msa_stats.get("local_alignment_path")
    alpha = curr_msa_stats["alpha"]
    n_branches = (2 * curr_msa_stats["n_seq"]) - 3
    if curr_msa_stats["use_spr_neighbours_training"]: #using original branch-lengths
        random_tree_per_site_ll_list, training_eval_running_time = raxml_compute_tree_per_site_ll(
            raxml_ll_eval_directory, local_file_path,
            random_trees_path, "sitelh_eval_brlen_opt", alpha=alpha, curr_msa_stats=curr_msa_stats,
            opt_brlen=False)
    elif brlen_generator_func is None:  #using optimized branch-lengths
        random_tree_per_site_ll_list, training_eval_running_time = raxml_compute_tree_per_site_ll(
            raxml_ll_eval_directory, local_file_path,
            random_trees_path, "sitelh_eval_brlen_opt", alpha=alpha, curr_msa_stats=curr_msa_stats,
            opt_brlen=True)
    else: #using other branch-lengths
        with open(random_trees_path, 'r') as RANDOM_TREES:
            random_trees_newick = RANDOM_TREES.read().split("\n")
        random_trees_objects = [generate_tree_object_from_newick(tree_newick) for tree_newick in random_trees_newick if
                                len(tree_newick) > 0]
        brlen_list_per_tree = [brlen_generator_func(n_branches, seed) for seed in
                               range(SEED, SEED + n_branches * len(random_trees_objects), n_branches)]
        random_trees_objects_with_brlen = [assign_brlen_to_tree_object(tree_obj, brlen_list) for tree_obj, brlen_list in
                                           zip(random_trees_objects, brlen_list_per_tree)]
        random_trees_newick_with_brlen = [tree_obj.write(format=1) for tree_obj in random_trees_objects_with_brlen]
        random_trees_with_brlen_path = os.path.join(curr_run_directory, "training_trees_with_brlen")
        with open(random_trees_with_brlen_path, 'w') as   RANDOM_TREES_WITH_BRLEN:
            RANDOM_TREES_WITH_BRLEN.write("\n".join(random_trees_newick_with_brlen))
        random_tree_per_site_ll_list, training_eval_running_time = raxml_compute_tree_per_site_ll(
            raxml_ll_eval_directory, local_file_path,
            random_trees_with_brlen_path, "sitelh_eval_w_brlen", alpha=alpha, curr_msa_stats=curr_msa_stats,
            opt_brlen=False)
    sitelh_df = pd.DataFrame(random_tree_per_site_ll_list, columns=list(range(len(random_tree_per_site_ll_list[0]))),
                             index=list(range(len(random_tree_per_site_ll_list))))
    sitelh_df.to_csv(output_csv_path, index=False)
    logging.debug(
        "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=output_csv_path))
    logging.info("Deleting dir content of {}".format(raxml_ll_eval_directory))
    # delete_dir_content(raxml_ll_eval_directory)
    return sitelh_df, training_eval_running_time
