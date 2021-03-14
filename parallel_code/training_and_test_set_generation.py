
from lasso_model_pipeline import *


def Lasso_training_and_test(brlen_generators, curr_msa_stats, training_size_options, random_trees_test_size):
    Lasso_folder = os.path.join(curr_msa_stats["curr_msa_version_folder"], "Lasso_folder")
    create_dir_if_not_exists(Lasso_folder)
    logging.info("Generating Lasso folder in {}".format(Lasso_folder))
    random_trees_folder = os.path.join(Lasso_folder, "random_tree_generation")
    curr_msa_stats["Lasso_folder"] = Lasso_folder
    curr_msa_stats["random_trees_folder"] = random_trees_folder
    create_dir_if_not_exists(random_trees_folder)
    random_trees_per_training_size = {}
    start_seed = SEED
    for training_size in training_size_options:
        logging.info("Generating common topologies for training size {} starting from seed {}".format(training_size,start_seed))
        training_random_tree_path_and_folder_list = generate_n_random_topologies_constant_brlen(training_size, random_trees_folder,
                                                                                                curr_msa_stats, "training", start_seed= start_seed)
        random_trees_per_training_size[training_size] = training_random_tree_path_and_folder_list
        start_seed += training_size
    logging.info(
        "Generating {} common topologies for common test set starting from seed {}".format(random_trees_test_size,start_seed))
    random_trees_test = generate_n_random_topologies_constant_brlen(random_trees_test_size, random_trees_folder, curr_msa_stats,
                                                     "test", start_seed=start_seed)
    test_folder = os.path.join(Lasso_folder, "test_{}_random_trees_eval".format(random_trees_test_size))
    create_dir_if_not_exists(test_folder)
    optimized_topologies_list_for_testing = generate_optimized_tree_topologies_for_testing(curr_msa_stats,random_trees_test ,test_folder)
    run_configurations = {}
    for brlen_generator_name in brlen_generators:
        brlen_run_directory = os.path.join(Lasso_folder, brlen_generator_name)
        create_dir_if_not_exists(brlen_run_directory)
        brlen_generator_func = brlen_generators.get(brlen_generator_name)
        brlen_start_seed = SEED
        logging.info(
            "Generating branch lengths for {} distribution starting from seed {}".format(brlen_generator_name, brlen_start_seed))
        for training_size in training_size_options:
            logging.info(
                "Generating branch lengths for training size {} starting from seed {}".format(training_size,
                                                                                             brlen_start_seed))
            training_size_directory = os.path.join(brlen_run_directory,
                                                   "training_{}_random_tree_eval".format(training_size))
            create_dir_if_not_exists(training_size_directory)
            training_sitelh, training_sitelh_path, last_seed_used = get_training_df(curr_msa_stats, brlen_generator_func,
                                                                    training_size_directory,
                                                                    random_trees_per_training_size[training_size],start_seed = brlen_start_seed)
            Lasso_results = apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats,
                                                                             curr_run_directory=training_size_directory,
                                                                             sitelh_training_df=training_sitelh,
                                                                             test_optimized_trees=optimized_topologies_list_for_testing)  # calculating positions_weight
            if brlen_generator_name not in run_configurations:
                run_configurations[brlen_generator_name] = {}
            run_configurations[brlen_generator_name][training_size] = Lasso_results
            brlen_start_seed = last_seed_used
    return run_configurations


def generate_n_random_topologies_constant_brlen(n, curr_run_directory, curr_msa_stats, name, start_seed):
    local_file_path = curr_msa_stats.get("local_alignment_path")
    random_tree_path_list= []
    basic_directory = os.path.join(curr_run_directory, "{name}_{n}".format(name=name,n=n))
    create_dir_if_not_exists(basic_directory)
    seed= start_seed
    for i in range(n):
        alpha = curr_msa_stats["alpha"]
        tree_folder = os.path.join(basic_directory, "random_tree_{i}".format(i=i))
        create_dir_if_not_exists(tree_folder)
        random_tree_generation_prefix = os.path.join(tree_folder, str(i))
        random_tree_path = generate_random_tree_topology_constant_brlen(alpha, local_file_path, random_tree_generation_prefix, seed=seed)
        random_tree_path_list.append(random_tree_path)
        seed=seed+1
    return random_tree_path_list




def generate_optimized_tree_topologies_for_testing(curr_msa_stats, random_trees_path_list,curr_run_directory):
    local_file_path = curr_msa_stats.get("local_alignment_path")
    optimized_random_trees_path_list = []
    for i, random_tree_path in enumerate(random_trees_path_list):
        prefix = "tree_{i}".format(i=i)
        optimized_tree = raxml_optimize_ll_on_given_tree_and_msa(local_file_path, prefix, random_tree_path, curr_msa_stats,
                                                                 curr_run_directory, weights=False, return_tree=True)
        optimized_random_trees_path_list.append(optimized_tree)
    return optimized_random_trees_path_list

def generate_per_site_ll_on_random_trees_for_training(curr_msa_stats, random_trees_path_list, brlen_generator_func, curr_run_directory, output_csv_path, start_seed):
    raxml_ll_eval_directory = os.path.join(curr_run_directory,"raxml_training_per_site_ll_eval")
    create_dir_if_not_exists(raxml_ll_eval_directory)
    local_file_path = curr_msa_stats.get("local_alignment_path")
    sitelh_ll_list = []
    alpha = curr_msa_stats["alpha"]
    seed = start_seed
    for i,random_tree_path in enumerate(random_trees_path_list):
        n_branches = (2*curr_msa_stats["n_seq"])-3
        prefix = "tree_{i}".format(i=i)
        tree_results_dir = os.path.join(raxml_ll_eval_directory,prefix)
        create_dir_if_not_exists(tree_results_dir)
        if brlen_generator_func is None:
            random_tree_per_site_ll_list = raxml_compute_tree_per_site_ll(tree_results_dir, local_file_path,
                                                                          random_tree_path, prefix, alpha,
                                                                          opt_brlen=True)
        else:
            tree_w_brlen_path = os.path.join(tree_results_dir,"tree_{i}_w_brlen.tree".format(i=i))
            assign_brlen(brlen_list=brlen_generator_func(size=n_branches, start_seed=seed),tree_path=random_tree_path,output_tree_path=tree_w_brlen_path)
            random_tree_per_site_ll_list = raxml_compute_tree_per_site_ll(tree_results_dir,  local_file_path,
                                                             tree_w_brlen_path, prefix, alpha,opt_brlen=False)
            seed=seed+n_branches
        sitelh_ll_list.append(random_tree_per_site_ll_list)
    sitelh_df = pd.DataFrame(sitelh_ll_list, columns=list(range(len(sitelh_ll_list[0]))),
                             index=list(range(len(sitelh_ll_list))))
    sitelh_df.to_csv(output_csv_path, index=False)
    logging.info(
        "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=output_csv_path))
    logging.info("Deleting dir content of {}".format(raxml_ll_eval_directory))
    delete_dir_content(raxml_ll_eval_directory)
    return sitelh_df, seed

#    delete_dir_content(random_trees_directory)





def get_training_df(curr_msa_stats, brlen_generator_func, curr_run_directory, random_trees_training, start_seed):
    training_output_csv_path = os.path.join(curr_run_directory,
                                            "training" + ".csv")
    logging.info("Generating training data in {}".format(training_output_csv_path))
    training_sitelh,last_seed_used = generate_per_site_ll_on_random_trees_for_training(curr_msa_stats=curr_msa_stats,
                                                                        random_trees_path_list=random_trees_training,
                                                                        brlen_generator_func=brlen_generator_func,
                                                                        curr_run_directory=curr_run_directory,
                                                                        output_csv_path=training_output_csv_path, start_seed = start_seed)
    return training_sitelh, training_output_csv_path, last_seed_used


def get_test_set_df(curr_msa_stats, brlen_generator_func, curr_run_directory, random_trees_test):
    test_output_csv_path = os.path.join(curr_run_directory,
                                        "test" + ".csv")
    logging.info("Generating test data in {}".format(test_output_csv_path))
    test_sitelh = generate_per_site_ll_on_random_trees_for_training(curr_msa_stats=curr_msa_stats,
                                                                    random_trees_path_list=random_trees_test,
                                                                    brlen_generator_func=brlen_generator_func,
                                                                    curr_run_directory=curr_run_directory,
                                                                    output_csv_path=test_output_csv_path)
    return test_sitelh, test_output_csv_path









