
from raxml import *
from spr_prune_and_regraft import *




def generate_n_random_topologies(n, curr_run_directory,curr_msa_stats,name):
    local_file_path = curr_msa_stats.get("local_alignment_path")
    random_tree_path_list= []
    basic_directory = os.path.join(curr_run_directory, "{name}_{n}".format(name=name,n=n))
    create_dir_if_not_exists(basic_directory)
    for i in range(n):
        alpha = curr_msa_stats["alpha"]
        tree_folder = os.path.join(basic_directory, "random_tree_{i}".format(i=i))
        create_dir_if_not_exists(tree_folder)
        random_tree_generation_prefix = os.path.join(tree_folder, str(i))
        random_tree_path=generate_random_tree_topology(alpha,local_file_path, random_tree_generation_prefix)
        random_tree_path_list.append(random_tree_path)
    return random_tree_path_list



def eval_per_site_ll_on_random_trees(curr_msa_stats, random_trees_path_list,brlen_generator_func, curr_run_directory, output_csv_path):
    raxml_ll_eval_directory = os.path.join(curr_run_directory,"raxml_tree_eval")
    create_dir_if_not_exists(raxml_ll_eval_directory)
    local_file_path = curr_msa_stats.get("local_alignment_path")
    n_trees = len(random_trees_path_list)
    logging.info("Generating " + str(n_trees) + " random trees ")
    sitelh_ll_list = []
    for i,random_tree_path in enumerate(random_trees_path_list):
        alpha = curr_msa_stats["alpha"]
        n_branches = 2*curr_msa_stats["n_seq"]-3
        prefix = "tree_{i}".format(i=i)
        if brlen_generator_func is None:
            create_dir_if_not_exists(prefix)
            random_tree_per_site_ll_list = raxml_compute_tree_per_site_ll(raxml_ll_eval_directory, local_file_path,
                                                                          random_tree_path, prefix, alpha,
                                                                          opt_brlen=True)
        else:
            tree_w_brlen_path = os.path.join(raxml_ll_eval_directory,"tree_{i}_w_brlen.tree".format(i=i))
            assign_brlen(brlen_list=brlen_generator_func(size=n_branches),tree_path=random_tree_path,output_tree_path=tree_w_brlen_path)
            random_tree_per_site_ll_list = raxml_compute_tree_per_site_ll(raxml_ll_eval_directory,  local_file_path,
                                                             tree_w_brlen_path, prefix, alpha,opt_brlen=False)
        sitelh_ll_list.append(random_tree_per_site_ll_list)
    sitelh_df = pd.DataFrame(sitelh_ll_list, columns=list(range(len(sitelh_ll_list[0]))),
                             index=list(range(len(sitelh_ll_list))))
    sitelh_df.to_csv(output_csv_path, index=False)
    logging.info(
        "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=output_csv_path))
    logging.info("Deleting dir content of {}".format(raxml_ll_eval_directory))
    delete_dir_content(raxml_ll_eval_directory)
    return sitelh_df

#    delete_dir_content(random_trees_directory)



def get_training_df(curr_msa_stats, brlen_generator_func, curr_run_directory, random_trees_training):
    training_output_csv_path = os.path.join(curr_run_directory,
                                            "training" + ".csv")
    logging.info("Generating training data in {}".format(training_output_csv_path))
    training_sitelh = eval_per_site_ll_on_random_trees(curr_msa_stats=curr_msa_stats,
                                                       random_trees_path_list=random_trees_training,
                                                       brlen_generator_func=brlen_generator_func,
                                                       curr_run_directory=curr_run_directory,
                                                       output_csv_path=training_output_csv_path)
    return training_sitelh, training_output_csv_path


def get_test_set_df(curr_msa_stats, brlen_generator_func, curr_run_directory, random_trees_test):
    test_output_csv_path = os.path.join(curr_run_directory,
                                        "test" + ".csv")
    logging.info("Generating test data in {}".format(test_output_csv_path))
    test_sitelh = eval_per_site_ll_on_random_trees(curr_msa_stats=curr_msa_stats,
                                                   random_trees_path_list=random_trees_test,
                                                   brlen_generator_func=brlen_generator_func,
                                                   curr_run_directory=curr_run_directory,
                                                   output_csv_path=test_output_csv_path)
    return test_sitelh, test_output_csv_path





