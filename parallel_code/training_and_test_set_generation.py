
from lasso_model_pipeline import *
import pickle


def Lasso_training_and_test(brlen_generators, curr_msa_stats, training_size_options, random_trees_test_size):
    Lasso_folder = os.path.join(curr_msa_stats["curr_msa_version_folder"], "Lasso_folder")
    create_dir_if_not_exists(Lasso_folder)
    logging.info("Generating Lasso folder in {}".format(Lasso_folder))
    random_trees_folder = os.path.join(Lasso_folder, "random_tree_generation")
    curr_msa_stats["Lasso_folder"] = Lasso_folder
    curr_msa_stats["random_trees_folder"] = random_trees_folder
    create_dir_if_not_exists(random_trees_folder)
    start_seed_random_trees = SEED
    if not curr_msa_stats["no_test_set"]:
        test_folder = os.path.join(Lasso_folder, "test_{}_random_trees_eval".format(random_trees_test_size))
        create_dir_if_not_exists(test_folder)
        logging.info("Generating test set based on {} random tree topologies".format(random_trees_test_size))
        test_random_trees_path, test_random_tree_generation_time = generate_n_random_topologies_constant_brlen(random_trees_test_size, random_trees_folder,
                                                                         curr_msa_stats,
                                                                         "test", seed= start_seed_random_trees)
        logging.info("Optimizing test set tree topologies".format(random_trees_test_size))
        optimized_test_topologies_path = generate_optimized_tree_topologies_for_testing(curr_msa_stats,
                                                                                        test_random_trees_path, test_folder)
        logging.info("Done with test set ")
    else:
        optimized_test_topologies_path = None

    run_configurations = {}
    for training_size in training_size_options:
        start_seed_random_trees += 1
        logging.info(
            "Generating common topologies for training size {} using seed {}".format(training_size, start_seed_random_trees))
        training_random_trees_path,training_tree_generation_elapsed_running_time = generate_n_random_topologies_constant_brlen(training_size, random_trees_folder,
                                                                                 curr_msa_stats, "training",
                                                                                 seed=start_seed_random_trees)
        logging.info("Done generating {} random topologies. It took {} seconds.".format(training_size,round(training_tree_generation_elapsed_running_time,2)))
        for brlen_generator_name in brlen_generators:
            logging.info("** Working on {} branch-length".format(brlen_generator_name))
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

            if os.path.exists(training_dump_baseline):
                logging.info("Using trainng results in {}".format(training_dump_baseline))
                with open(training_dump_baseline, 'rb') as handle:
                    training_results = pickle.load(handle)
                    training_sitelh, training_eval_time = training_results["training_sitelh"], training_results[
                        "training_eval_time"]
            training_sitelh, training_sitelh_path,training_eval_time = get_training_df(curr_msa_stats, brlen_generator_func,
                                                                    training_size_directory,
                                                                                    training_random_trees_path)
            logging.info("Done evaluating topologies based on {} branch lengths. It took {} seconds".format(brlen_generator_name,training_eval_time))

            logging.info('Applying Lasso on current training data')
            Lasso_results = apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats,
                                                                             curr_run_directory=training_size_directory,
                                                                             sitelh_training_df=training_sitelh,
                                                                             test_optimized_trees_path=optimized_test_topologies_path)  # calculating positions_weight
            Lasso_results.update({'training_random_trees_generation_time':training_tree_generation_elapsed_running_time ,
                                  'training_evaluation_time' : training_eval_time
                                  })
            logging.info("Lasso results: \n {}".format({k:  Lasso_results[k] for k in  Lasso_results.keys() if
                                                          k not in ["lasso_training_X","lasso_training_Y","lasso_chosen_locis","lasso_coeffs","lasso_chosen_weights"]
                                                          }))
            if brlen_generator_name not in run_configurations:
                run_configurations[brlen_generator_name] = {}
            run_configurations[brlen_generator_name][training_size] = Lasso_results
    return run_configurations


def generate_n_random_topologies_constant_brlen(n, curr_run_directory, curr_msa_stats, name, seed):
    local_file_path = curr_msa_stats.get("local_alignment_path")
    basic_directory = os.path.join(curr_run_directory, "{name}_{n}".format(name=name,n=n))
    create_dir_if_not_exists(basic_directory)
    seed= seed
    alpha = curr_msa_stats["alpha"]
    rrandom_tree_path,elapsed_running_time= generate_n_random_tree_topology_constant_brlen(n,alpha, local_file_path, basic_directory,curr_msa_stats, seed=seed)
    return rrandom_tree_path,elapsed_running_time




def generate_optimized_tree_topologies_for_testing(curr_msa_stats, test_random_trees_path,curr_run_directory):
    local_file_path = curr_msa_stats.get("local_alignment_path")
    baseline_trees_path = os.path.join(curr_run_directory, "test_opt","test_opt.raxml.mlTrees").replace(curr_msa_stats["run_prefix"],
                                                             curr_msa_stats["training_set_baseline_run_prefix"])
    if os.path.exists(baseline_trees_path):
        logging.info("Using test set in {training_set_baseline_run_prefix} ; path is: {baseline_trees_path}")
        return baseline_trees_path
    else:
        logging.info("Generating test set from beggining")
    optimized_trees_path = raxml_optimize_trees_for_given_msa(local_file_path, "test_opt", test_random_trees_path, curr_msa_stats,
                                                              curr_run_directory, weights=False, return_trees_file=True)

    return optimized_trees_path

def generate_per_site_ll_on_random_trees_for_training(curr_msa_stats, random_trees_path, brlen_generator_func, curr_run_directory, output_csv_path):
    raxml_ll_eval_directory = os.path.join(curr_run_directory,"raxml_training_per_site_ll_eval")
    create_dir_if_not_exists(raxml_ll_eval_directory)
    local_file_path = curr_msa_stats.get("local_alignment_path")
    alpha = curr_msa_stats["alpha"]
    n_branches = (2*curr_msa_stats["n_seq"])-3
    if brlen_generator_func is None:
        random_tree_per_site_ll_list, training_eval_running_time = raxml_compute_tree_per_site_ll(raxml_ll_eval_directory, local_file_path,
                                                                      random_trees_path, "sitelh_eval_brlen_opt", alpha=alpha,curr_msa_stats=curr_msa_stats,
                                                                      opt_brlen=True)
    else:
        with open(random_trees_path,'r') as RANDOM_TREES:
            random_trees_newick = RANDOM_TREES.read().split("\n")
        random_trees_objects = [generate_tree_object_from_newick(tree_newick) for tree_newick in random_trees_newick if len(tree_newick)>0]
        brlen_list_per_tree = [brlen_generator_func(n_branches,seed) for seed in range(SEED,SEED+n_branches*len(random_trees_objects),n_branches)]
        random_trees_objects_with_brlen = [assign_brlen_to_tree_object(tree_obj,brlen_list) for tree_obj,brlen_list in zip(random_trees_objects,brlen_list_per_tree)]
        random_trees_newick_with_brlen = [tree_obj.write(format=1) for tree_obj in   random_trees_objects_with_brlen]
        random_trees_with_brlen_path = os.path.join(curr_run_directory,"training_trees_with_brlen")
        with open( random_trees_with_brlen_path,'w') as   RANDOM_TREES_WITH_BRLEN:
            RANDOM_TREES_WITH_BRLEN.write("\n".join(random_trees_newick_with_brlen))
        random_tree_per_site_ll_list,training_eval_running_time = raxml_compute_tree_per_site_ll(raxml_ll_eval_directory,  local_file_path,
                                                                      random_trees_with_brlen_path, "sitelh_eval_w_brlen", alpha=alpha,curr_msa_stats=curr_msa_stats,opt_brlen=False)
    sitelh_df = pd.DataFrame(random_tree_per_site_ll_list, columns=list(range(len(random_tree_per_site_ll_list[0]))),
                             index=list(range(len(random_tree_per_site_ll_list))))
    sitelh_df.to_csv(output_csv_path, index=False)
    logging.info(
        "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=output_csv_path))
    logging.info("Deleting dir content of {}".format(raxml_ll_eval_directory))
    #delete_dir_content(raxml_ll_eval_directory)
    return sitelh_df,training_eval_running_time





def get_training_df(curr_msa_stats, brlen_generator_func, curr_run_directory, random_trees_training):
    training_output_csv_path = os.path.join(curr_run_directory,
                                            "training" + ".csv")
    training_dump = os.path.join(curr_run_directory, 'training_set.dump')
    training_dump_baseline = training_dump.replace(curr_msa_stats["run_prefix"],
                                                                     curr_msa_stats["training_set_baseline_run_prefix"])

    if os.path.exists(training_dump_baseline):
        logging.info("Using trainng results in {}".format(training_dump_baseline))
        with open(training_dump_baseline, 'rb') as handle:
            training_results = pickle.load(handle)
            training_sitelh,training_eval_time = training_results["training_sitelh"], training_results["training_eval_time"]
    else:
        training_sitelh, training_eval_time = generate_per_site_ll_on_random_trees_for_training(curr_msa_stats=curr_msa_stats,
                                                                                               random_trees_path=random_trees_training,
                                                                                               brlen_generator_func=brlen_generator_func,
                                                                                               curr_run_directory=curr_run_directory,
                                                                                               output_csv_path=training_output_csv_path)
        training_results = {"training_sitelh" : training_sitelh, "training_eval_time": training_eval_time}
        with open(training_dump, 'wb') as handle:
            pickle.dump(training_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return training_sitelh, training_output_csv_path,training_eval_time












