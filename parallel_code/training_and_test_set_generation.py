
from raxml import *
from spr_prune_and_regraft import *




def generate_site_lh_data(curr_msa_stats, n_iter,brlen_generator_func,curr_run_directory,output_csv_path):
    create_dir_if_not_exists(curr_run_directory)
    local_file_path = curr_msa_stats.get("local_alignment_path")
    logging.info("Generating " + str(n_iter) + " random trees ")
    sitelh_ll_list = []
    random_trees_directory = os.path.join(curr_run_directory,"sitelh_generation")
    create_dir_if_not_exists(random_trees_directory)
    for i in range(n_iter):
        alpha = curr_msa_stats["alpha"]
        random_tree_generation_prefix = os.path.join(random_trees_directory, str(i))
        random_tree_path = generate_random_tree_topology(alpha,  local_file_path, random_tree_generation_prefix)
        n_branches = 2*curr_msa_stats["n_seq"]-3
        if brlen_generator_func is None:
            random_tree_per_site_ll_list = raxml_compute_tree_per_site_ll(random_trees_directory, local_file_path,
                                                                          random_tree_path, str(i), alpha,
                                                                          opt_brlen=True)
        else:
            assign_brlen(brlen_list=brlen_generator_func(size=n_branches),tree_path=random_tree_path)
            random_tree_per_site_ll_list = raxml_compute_tree_per_site_ll(random_trees_directory,  local_file_path,
                                                             random_tree_path, str(i), alpha,opt_brlen=False)
        sitelh_ll_list.append(random_tree_per_site_ll_list)
    sitelh_df = pd.DataFrame(sitelh_ll_list, columns=list(range(len(sitelh_ll_list[0]))),
                             index=list(range(len(sitelh_ll_list))))
    sitelh_df.to_csv(output_csv_path, index=False)
    delete_dir_content(random_trees_directory)
    logging.info(
        "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=output_csv_path))
    return sitelh_df









