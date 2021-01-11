
from raxml import *
from spr_prune_and_regraft import *




def generate_site_lh_data(curr_msa_stats, n_iter,name,brlen_generator_func,curr_run_directory):
    create_dir_if_not_exists(curr_run_directory)
    csv_path = os.path.join(curr_run_directory, curr_msa_stats.get("file_name")+"_"+name +".csv")
    csv_path_baseline = csv_path.replace(curr_msa_stats["run_prefix"],curr_msa_stats["baseline_run_prefix"])
    if os.path.exists(csv_path_baseline):
        baseline_df = pd.read_csv(csv_path_baseline)
        if not baseline_df.empty:
            logging.info("Using baseline DataFrame in {csv_path_baseline}, and copying it to {csv_path}".foramt(csv_path_baseline=csv_path_baseline, csv_path=csv_path))
            sitelh_df = baseline_df
            sitelh_df.to_csv(csv_path, index=False)
            return sitelh_df
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
    sitelh_df.to_csv(csv_path, index=False)
    delete_dir_content(random_trees_directory)
    logging.info(
        "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=csv_path))
    return sitelh_df









