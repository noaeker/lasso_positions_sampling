from raxml import *
import numpy as np
from scipy import stats
from spr_prune_and_regraft import *
from sklearn.metrics import *


def write_spr_log_message(spr_log_file_object, rgrft_path, best_ll, ll, best_topology_path):
    if ll > best_ll:
        spr_log_file_object.write("Found a better tree!!!! in " + rgrft_path + "\n")
        spr_log_file_object.write("    1. best ll is: " + str(best_ll) + "\n")
        spr_log_file_object.write("    2. current tree in " + rgrft_path + " has log likelihood of " + str(ll) + "\n")
        spr_log_file_object.write("    copying topology in " + rgrft_path + " to " + best_topology_path + "\n")
        spr_log_file_object.write("updating best likelihood to be " + str(ll) + "\n")
    elif ll > -np.infty:
        spr_log_file_object.write("Not found a better tree in " + rgrft_path + "\n")
        spr_log_file_object.write("    1. best ll is: " + str(best_ll) + "\n")
        spr_log_file_object.write("    2. current tree in " + rgrft_path + " has log likelihood of " + str(ll) + "\n")


def get_true_ll_values_and_sitelh(weights_file_path, curr_msa_stats, trees_path, curr_run_directory, trees_ll):
    if weights_file_path:
        if not curr_msa_stats["compute_all_true_ll"]:
            return [-1] * len(trees_ll), pd.DataFrame()
        if curr_msa_stats["compute_per_site_ll_values"]:
            random_tree_per_site_ll_list, random_tree_per_site_ll_list_eval_time = raxml_compute_tree_per_site_ll(
                curr_run_directory, curr_msa_stats["local_alignment_path"],
                trees_path, "rgrft_ll_eval_on_full_MSA", alpha=curr_msa_stats["alpha"],
                curr_msa_stats=curr_msa_stats,
                opt_brlen=True)
            true_sitelh_df = pd.DataFrame(random_tree_per_site_ll_list,
                                          columns=list(range(len(random_tree_per_site_ll_list[0]))),
                                          index=list(range(len(random_tree_per_site_ll_list))))
            trees_true_ll = list(true_sitelh_df.sum(axis=1))
        else:
            trees_true_ll, trees_true_optimized_objects, time_rgft_eval_true = raxml_optimize_trees_for_given_msa(
                curr_msa_stats["local_alignment_path"],
                "rgrft_ll_eval_on_full_MSA", trees_path,
                curr_msa_stats, curr_run_directory,
                weights=None, n_cpus=curr_msa_stats["n_cpus_full"])
        return trees_true_ll, pd.DataFrame()
    else:
        return trees_ll, pd.DataFrame()


def compute_true_ll_of_best_tree_of_spr_iteration(weights_file_path, best_tree_object, curr_run_directory, curr_msa_stats,
                                                  trees_true_ll, best_ll_index, best_ll, top_x_to_test):
    if weights_file_path and top_x_to_test == 1:
        if curr_msa_stats["compute_all_true_ll"]:
            return trees_true_ll[best_ll_index]
        else:
            best_first_phase_newick = (best_tree_object.write(format=1))
            best_tree_path = os.path.join(curr_run_directory, "iteration_best_spr_tree")
            with open(best_tree_path, 'w') as BEST_TREE:
                BEST_TREE.write(best_first_phase_newick)
            best_tree_true_ll, best_tree_true_optimized_object, best_tree_true_eval_true = raxml_optimize_trees_for_given_msa(
                curr_msa_stats["local_alignment_path"],
                "best_iter_tree_eval_full_MSA", best_tree_path,
                curr_msa_stats, curr_run_directory,
                weights=None, n_cpus=curr_msa_stats["n_cpus_full"])
            best_true_ll = best_tree_true_ll
            return best_true_ll
    else:
        return best_ll


def regression_correct_lasso_ll_values(lasso_intercept,weights_file_path,trees_ll):
    if weights_file_path:
        ll_fixed = [((ll) + lasso_intercept) / INTEGER_CONST for ll in trees_ll]
        return ll_fixed
    else:
        return trees_ll




def write_tree_objects_to_file(trees_objects,curr_run_directory, top_trees_file_name):
    top_ll_tree_objects = np.array(trees_objects)
    top_ll_trees_newick = "\n".join([tree_object.write(format=1) for tree_object in top_ll_tree_objects])
    top_ll_trees_path = os.path.join(curr_run_directory, top_trees_file_name)
    with open(top_ll_trees_path, 'w') as TOP_LL_TREES:
        TOP_LL_TREES.write(top_ll_trees_newick)
    return top_ll_trees_path





def get_non_greedy_optimized_SPR_neighbours(curr_msa_stats, MSA_path, unique_trees_path, weights_file_path,lasso_intercept, curr_run_directory, n_cpus):
    if curr_msa_stats["optimized_neighbours_per_iter"] > 1:
        logging.info("Evaluating (no brlen opt) LL of all SPR neighbours")
        trees_ll_no_brlen, trees_optimized_objects_no_brlen, time_rgft_eval_no_brlen = raxml_optimize_trees_for_given_msa(
            MSA_path, "rgrft_ll_eval_no_brlen",
            unique_trees_path,
            curr_msa_stats,
            curr_run_directory, opt_brlen=False,
            weights=weights_file_path, n_cpus= n_cpus)
        indices_of_spr_candidates_for_brlen_opt = (-np.array(trees_ll_no_brlen)).argsort()[
                                                  :curr_msa_stats["optimized_neighbours_per_iter"]]
        tree_objects_of_spr_candidates_for_brlen_opt = np.array(trees_optimized_objects_no_brlen)[
            indices_of_spr_candidates_for_brlen_opt]
        spr_candidates_for_brlen_opt_file = write_tree_objects_to_file(tree_objects_of_spr_candidates_for_brlen_opt,
                                                                       curr_run_directory,
                                                                       "spr_candidates_for_brlen_opt.trees")
        logging.info("About to optimize LL of most promising {t} topologies".format(t=curr_msa_stats["optimized_neighbours_per_iter"]))
    else:
        logging.info("Fully optimizing all SPR neighbours")
        spr_candidates_for_brlen_opt_file = unique_trees_path
    ll_spr_candidates_for_brlen, optimized_objects_spr_candidates_for_brlen, time_rgft_eval_true_spr_candidates_for_brlen = raxml_optimize_trees_for_given_msa(
        MSA_path, "rgrft_ll_eval_brlen",
        spr_candidates_for_brlen_opt_file,
        curr_msa_stats,
        curr_run_directory,
        weights=weights_file_path ,n_cpus=n_cpus)
    ll_spr_candidates_for_brlen_corrected = regression_correct_lasso_ll_values(lasso_intercept,weights_file_path,
                                                                               ll_spr_candidates_for_brlen)
    return ll_spr_candidates_for_brlen_corrected, optimized_objects_spr_candidates_for_brlen, time_rgft_eval_true_spr_candidates_for_brlen,time_rgft_eval_no_brlen


def re_optimize_some_SPR_neighbours_no_weights(ll_spr_candidates_for_brlen_corrected,top_x_true_trees,optimized_objects_spr_candidates_for_brlen, curr_msa_stats, curr_run_directory):
    top_ll_indices = (-np.array(ll_spr_candidates_for_brlen_corrected)).argsort()[:top_x_true_trees]
    top_ll_tree_objects = np.array(optimized_objects_spr_candidates_for_brlen)[top_ll_indices]
    top_ll_trees_newick = "\n".join([tree_object.write(format=1) for tree_object in top_ll_tree_objects])
    top_ll_trees_path = os.path.join(curr_run_directory, "lasso_top_ll_trees_file.trees")
    with open(top_ll_trees_path, 'w') as TOP_LL_TREES:
        TOP_LL_TREES.write(top_ll_trees_newick)
    top_trees_true_ll, top_trees_true_optimized_objects, time_rgft_eval_true = raxml_optimize_trees_for_given_msa(
        curr_msa_stats["local_alignment_path"],
        "lasso_re_optimization_on_full_MSA", top_ll_trees_path,
        curr_msa_stats, curr_run_directory,
        weights=None, n_cpus= curr_msa_stats["n_cpus_full"])  # optimize without weights
    best_ll = max(top_trees_true_ll)
    best_ll_index = top_trees_true_ll.index(best_ll)
    best_tree_object = top_trees_true_optimized_objects[best_ll_index]
    return best_ll, best_ll_index, best_tree_object,time_rgft_eval_true


# def SPR_greedy_iteration(curr_best_ll, MSA_path, curr_msa_stats, starting_tree_object,
#                   curr_run_directory,
#                   lasso_weights = None):
#     starting_tree_spr_neighbours = get_possible_spr_moves(starting_tree_object, rearr_dist=curr_msa_stats["rearr_dist"])
#     regrafted_trees = [generate_neighbour(starting_tree_object, spr_neighbour) for spr_neighbour in
#                        starting_tree_spr_neighbours]
#     regrafted_trees_newick =[regrafted_tree.write(format=1) for regrafted_tree in regrafted_trees]
#     trees_path = os.path.join(curr_run_directory, "iteration_spr_trees")
#     time_rgft_eval_no_brlen = 0
#     time_rgft_eval_brlen = 0
#     no_brlen_data = {}
#     n_no_brlen_eval = 0
#     n_brlen_eval =0
#     for tree in regrafted_trees_newick:
#         with open(trees_path, 'w') as TREES:
#             TREES.write(tree)
#         n_no_brlen_eval+=1
#         tree_ll_no_brlen, tree_optimized_objects_no_brlen, tree_time_rgft_eval_no_brlen = raxml_optimize_trees_for_given_msa(
#             MSA_path, "rgrft_ll_eval_no_brlen",
#             trees_path,
#             curr_msa_stats,
#             curr_run_directory, opt_brlen=False,
#             weights=lasso_weights)
#         time_rgft_eval_no_brlen = time_rgft_eval_no_brlen + tree_time_rgft_eval_no_brlen
#         no_brlen_data[tree_ll_no_brlen] = tree_optimized_objects_no_brlen
#         if tree_ll_no_brlen> curr_best_ll:
#             return tree_optimized_objects_no_brlen, tree_ll_no_brlen, pd.DataFrame(), pd.DataFrame(), time_rgft_eval_no_brlen, time_rgft_eval_no_brlen, time_rgft_eval_brlen, 0, n_no_brlen_eval, n_brlen_eval
#     else:
#         max_neighbours_with_brlen = curr_msa_stats["optimized_neighbours_per_iter"]
#         optimized_tree_path = os.path.join(curr_run_directory, "iteration_spr_trees")
#         for ll, tree_object in sorted(no_brlen_data.items(),reverse=True):
#             with open( optimized_tree_path, 'w') as TREES:
#                 TREES.write(tree_object.write(format=1))
#                 n_brlen_eval += 1
#             tree_ll, tree_optimized_object, tree_time_rgft_eval = raxml_optimize_trees_for_given_msa(
#                 MSA_path, "rgrft_ll_eval_no_brlen",
#                 trees_path,
#                 curr_msa_stats,
#                 curr_run_directory, opt_brlen=True,
#                 weights=lasso_weights)
#             time_rgft_eval_brlen = time_rgft_eval_brlen+tree_time_rgft_eval
#             if max_neighbours_with_brlen==0 or tree_ll> curr_best_ll:
#                 overall_iteration_time = time_rgft_eval_brlen + time_rgft_eval_no_brlen
#                 best_tree_true_ll, best_tree_true_optimized_object, best_tree_true_eval_true = raxml_optimize_trees_for_given_msa(
#                 curr_msa_stats["local_alignment_path"],
#                 "best_iter_tree_eval_full_MSA", best_tree_path,
#                 curr_msa_stats, curr_run_directory,
#                 weights=None)
#                 results_dict = {"best_tree_object": tree_optimized_object, "best_ll": tree_ll, "best_true_ll": best_true_ll,
#                                 "ll_comparison_df": pd.DataFrame(), "true_sitelh_df": pd.DataFrame(),
#                                 "iteration_time": overall_iteration_time, "iteration_time_brlen": time_rgft_eval_brlen,
#                                 "iteration_time_no_brlen": time_rgft_eval_no_brlen,
#                                 "re_optimization_time": 0}
#                 return results_dict
#                 return tree_optimized_object,tree_ll,pd.DataFrame(), pd.DataFrame(),overall_iteration_time,time_rgft_eval_no_brlen, time_rgft_eval_brlen, 0 , n_no_brlen_eval, n_brlen_eval
#             max_neighbours_with_brlen = max_neighbours_with_brlen -1



    #return best_tree_object, best_ll, best_true_ll, ll_comparison_df, true_sitelh_df, iteration_time,time_rgft_eval_true_spr_candidates_for_brlen,time_rgft_eval_no_brlen,re_optimization_time

def SPR_iteration(iteration_number, MSA_path, curr_msa_stats, starting_tree_object,
                  curr_run_directory,
                  weights_file_path,lasso_intercept, top_x_true_trees, n_cpus):
    iteration_time = 0
    re_optimization_time=0
    add_internal_names(starting_tree_object)
    starting_tree_object.get_tree_root().name = "ROOT"
    logging.debug(str(starting_tree_object.write(format=1)) + "\n")
    starting_tree_spr_neighbours = get_possible_spr_moves(starting_tree_object, rearr_dist = curr_msa_stats["rearr_dist"])
    regrafted_trees = [generate_neighbour(starting_tree_object, spr_neighbour) for spr_neighbour in
                       starting_tree_spr_neighbours]
    regrafted_trees_newick = "\n".join([regrafted_tree.write(format=1) for regrafted_tree in regrafted_trees])
    trees_path = os.path.join(curr_run_directory, "iteration_spr_trees")
    with open(trees_path, 'w') as TREES:
        TREES.write(regrafted_trees_newick)
    unique_trees_path = filter_unique_topologies(curr_run_directory, trees_path, len(regrafted_trees))
    ll_spr_candidates_for_brlen_corrected, optimized_objects_spr_candidates_for_brlen, iteration_time_brlen,iteration_time_no_brlen = get_non_greedy_optimized_SPR_neighbours(curr_msa_stats, MSA_path, unique_trees_path, weights_file_path,lasso_intercept, curr_run_directory, n_cpus)
    iteration_time+= iteration_time_brlen+iteration_time_no_brlen
    if top_x_true_trees > 1 and weights_file_path:
        logging.debug(f"SPR iteration {iteration_number} : testing {top_x_true_trees} best Lasso tree objects")
        best_ll, best_ll_index, best_tree_object, re_optimization_time = re_optimize_some_SPR_neighbours_no_weights(ll_spr_candidates_for_brlen_corrected,top_x_true_trees,optimized_objects_spr_candidates_for_brlen, curr_msa_stats, curr_run_directory)
        iteration_time += re_optimization_time
    else:
        best_ll = max(ll_spr_candidates_for_brlen_corrected)
        best_ll_index = ll_spr_candidates_for_brlen_corrected.index(best_ll)
        best_tree_object = optimized_objects_spr_candidates_for_brlen[best_ll_index]
    trees_true_ll, true_sitelh_df = get_true_ll_values_and_sitelh(weights_file_path, curr_msa_stats, trees_path,
                                                                  curr_run_directory, ll_spr_candidates_for_brlen_corrected)
    best_true_ll = compute_true_ll_of_best_tree_of_spr_iteration(weights_file_path, best_tree_object, curr_run_directory,
                                                                 curr_msa_stats, trees_true_ll, best_ll_index, best_ll,
                                                                 top_x_true_trees)
    ll_comparison_df = pd.DataFrame(
        {'full msa ll': trees_true_ll, 'sampled msa ll': ll_spr_candidates_for_brlen_corrected, 'iteration number': iteration_number}
    )

    results_dict = {"best_tree_object":best_tree_object, "best_ll" : best_ll,"best_true_ll": best_true_ll, "ll_comparison_df":ll_comparison_df,"true_sitelh_df": true_sitelh_df, "iteration_time": iteration_time, "iteration_time_brlen" : iteration_time_brlen, "iteration_time_no_brlen":iteration_time_no_brlen, "re_optimization_time":re_optimization_time, "n_neighbours": len(regrafted_trees) }
    return results_dict



def get_true_and_sampled_starting_tree_ll(reduced_MSA_path, run_unique_name, starting_tree_path, curr_msa_stats,
                                          curr_run_directory, weights_file_path,lasso_intercept, n_cpus):
    search_starting_tree_ll, tree_objects, elapsed_running_time_starting_eval = raxml_optimize_trees_for_given_msa(
        reduced_MSA_path,
        "starting_tree_ll_eval_" + run_unique_name,
        starting_tree_path,
        curr_msa_stats,
        curr_run_directory=curr_run_directory,
        weights=weights_file_path, n_cpus = n_cpus)
    if weights_file_path:
        search_starting_tree_ll = regression_correct_lasso_ll_values(lasso_intercept,weights_file_path, [search_starting_tree_ll])[0]
        search_true_starting_tree_ll, tree_objects, elapsed_running_time_starting_eval_true = raxml_optimize_trees_for_given_msa(
            curr_msa_stats["local_alignment_path"], "starting_tree_ll_eval_on_full_" + run_unique_name,
            starting_tree_path,
            curr_msa_stats, curr_run_directory=curr_run_directory, weights=None, n_cpus= curr_msa_stats["n_cpus_Lasso"])
    else:
        search_true_starting_tree_ll = search_starting_tree_ll

    return search_starting_tree_ll, search_true_starting_tree_ll


def SPR_search(MSA_path, run_unique_name, curr_msa_stats, starting_tree_path, starting_tree_object,
               curr_run_directory,
               weights_file_path,lasso_intercept, top_x_true_trees, starting_tree_ll=None, n_cpus = 1):
    ll_comparison_df = pd.DataFrame()
    running_times_per_iter = []
    no_brlen_times_per_iter = []
    brlen_per_iter = []
    re_optimization_time_per_iter = []
    actual_search_training_df = pd.DataFrame()
    spr_iterations_performed_so_far = 0
    total_spr_neighbours_evaluated=0
    search_starting_tree_ll, search_true_starting_tree_ll = get_true_and_sampled_starting_tree_ll(MSA_path,
                                                                                                  run_unique_name,
                                                                                                  starting_tree_path,
                                                                                                  curr_msa_stats,
                                                                                                  curr_run_directory,
                                                                                                  weights_file_path,lasso_intercept, n_cpus= n_cpus)

    if top_x_true_trees > 1:
        search_starting_tree_ll = search_true_starting_tree_ll
    LL_per_iteration_list = [search_starting_tree_ll]
    TRUE_LL_per_iteration_list = [search_true_starting_tree_ll]
    logging.debug("Search starting tree ll = {} Search true starting tree ll = {}".format(search_starting_tree_ll,
                                                                                          search_true_starting_tree_ll))
    curr_best_tree_ll, curr_best_tree_true_ll = search_starting_tree_ll, search_true_starting_tree_ll
    if starting_tree_object is None:
        starting_tree_object = generate_tree_object_from_newick(starting_tree_path)
    curr_best_tree_object = starting_tree_object
    while True:
        curr_iter_run_directory = os.path.join(curr_run_directory, "iter_" + str(spr_iterations_performed_so_far))
        create_or_clean_dir(curr_iter_run_directory)
        logging.debug("iteration number: " + str(spr_iterations_performed_so_far))
        new_iteration_results= SPR_iteration(
            spr_iterations_performed_so_far, MSA_path, curr_msa_stats, curr_best_tree_object,
            curr_iter_run_directory,
            weights_file_path, lasso_intercept,
            top_x_true_trees, n_cpus
        )
        logging.debug(
            "Our current best tree ll is {} (its true ll is {}), best neighbour ll is {}".format(curr_best_tree_ll,
                                                                                                 curr_best_tree_true_ll,
                                                                                                 new_iteration_results["best_ll"]))
        ll_comparison_df = ll_comparison_df.append(new_iteration_results["ll_comparison_df"])
        running_times_per_iter.append(new_iteration_results["iteration_time"])
        brlen_per_iter.append(new_iteration_results["iteration_time_brlen"])
        no_brlen_times_per_iter.append(new_iteration_results["iteration_time_no_brlen"])
        re_optimization_time_per_iter.append(new_iteration_results["re_optimization_time"])

        actual_search_training_df = actual_search_training_df.append(new_iteration_results["true_sitelh_df"])
        if new_iteration_results["best_ll"] - curr_best_tree_ll <= EPSILON:
            logging.debug(
                "Difference between best spr neighbour and current tree <= {}, stopping SPR search\n".format(EPSILON))
            break
        logging.debug("Updating best neighbour to be our current best tree! ")
        ### Updating current iteration results and preparing for next iteration:
        spr_iterations_performed_so_far = spr_iterations_performed_so_far + 1
        curr_best_tree_object = new_iteration_results["best_tree_object"]
        curr_best_tree_ll = new_iteration_results["best_ll"]
        curr_best_tree_true_ll = new_iteration_results["best_true_ll"]
        LL_per_iteration_list += [curr_best_tree_ll]
        TRUE_LL_per_iteration_list += [curr_best_tree_true_ll]
        total_spr_neighbours_evaluated = total_spr_neighbours_evaluated+new_iteration_results["n_neighbours"]
    curr_best_tree_path = os.path.join(curr_run_directory, "search_best_tree_path")
    with open(curr_best_tree_path,'w') as BEST_TREE:
        BEST_TREE.write(curr_best_tree_object.write(format=1))
    search_results = {
        "search_best_ll": curr_best_tree_ll,
        "search_starting_tree_ll": starting_tree_ll,
        "search_best_true_ll": curr_best_tree_true_ll,
        "search_best_topology_newick": curr_best_tree_object.write(format=1),
        "search_starting_tree_newick": starting_tree_object.write(format=1),
        "ll_comparison_df": ll_comparison_df,
        "actual_search_training_df":  actual_search_training_df,
        "ll_per_iteration_list": LL_per_iteration_list,
        "TRUE_ll_per_iteration_list": TRUE_LL_per_iteration_list,
        "search_best_tree_object": curr_best_tree_object,
        "search_best_tree_path": curr_best_tree_path,
        "search_spr_moves": spr_iterations_performed_so_far,
        "running_time_per_iter": running_times_per_iter,
        "total_search_running_time": sum(running_times_per_iter),
        "total_brlen_time" : sum(brlen_per_iter) ,
        "total_no_brlen_time" : sum(no_brlen_times_per_iter),
        "total_reoptimization_time" : sum (re_optimization_time_per_iter),
        "total_spr_neighbours_evaluated": total_spr_neighbours_evaluated
    }

    return search_results

def analyze_ll_comparison_df(ll_comparison_df):
    mistake_cnt = 0
    for iteration in ll_comparison_df['iteration number'].unique():
        curr_iter_ll_comparison_df = ll_comparison_df[ll_comparison_df['iteration number'] == iteration]
        maxvalueIndexLabel = curr_iter_ll_comparison_df.idxmax()
        if maxvalueIndexLabel['full msa ll'] != maxvalueIndexLabel['sampled msa ll']:
            mistake_cnt += 1
    try:
        rho_pearson, pval_pearson = stats.pearsonr(ll_comparison_df['full msa ll'],
                                                   ll_comparison_df['sampled msa ll'])
        rho_spearman, pval_spearman = stats.spearmanr(ll_comparison_df['full msa ll'],
                                                      ll_comparison_df['sampled msa ll'])
        mse = mean_squared_error(ll_comparison_df['full msa ll'], ll_comparison_df['sampled msa ll'])
    except:
        rho_pearson, pval_pearson, mistake_cnt, rho_spearman, pval_spearman, mse = -1, -1, -1, -1

    return rho_pearson, pval_pearson, rho_spearman, pval_spearman, mse, mistake_cnt


def SPR_analysis(lasso_configurations,SPR_chosen_starting_tree_path, curr_msa_stats, curr_run_directory,
                 full_run=False):
    run_unique_name = "spr"
    if not os.path.exists(curr_run_directory):
        os.mkdir(curr_run_directory)
    if (full_run):
        full_data_param_dict = SPR_search(
            MSA_path=curr_msa_stats["local_alignment_path"],
            run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=SPR_chosen_starting_tree_path,
            starting_tree_object=None,
            curr_run_directory=curr_run_directory, weights_file_path=None,lasso_intercept= -1, top_x_true_trees=1, n_cpus= curr_msa_stats["n_cpus_full"])
        full_data_SPR_result = {"naive_SPR_ll": full_data_param_dict["search_best_ll"],
                                "naive_SPR_starting_tree" : SPR_chosen_starting_tree_path,
                                "naive_SPR_spr_moves": full_data_param_dict["search_spr_moves"],
                                "naive_SPR_tree_newick": full_data_param_dict["search_best_topology_newick"]
            , "naive_SPR_ll_per_iteration": full_data_param_dict["ll_per_iteration_list"],
                                "SPR_search_starting_tree_ll": full_data_param_dict["search_starting_tree_ll"],
                                "SPR_search_starting_tree_newick": full_data_param_dict["search_starting_tree_newick"],
                                "naive_spr_brlen_running_time" : full_data_param_dict["total_brlen_time"],
                                "naive_spr_no_brlen_running_time": full_data_param_dict["total_no_brlen_time"],
                                "naive_SPR_running_time": full_data_param_dict["total_search_running_time"],
                                "naive_SPR_total_spr_neighbours_evaluated": full_data_param_dict["total_spr_neighbours_evaluated"]


                                }
        full_iterations_data = pd.DataFrame({"ll": full_data_SPR_result["naive_SPR_ll_per_iteration"],"true_ll": full_data_SPR_result["naive_SPR_ll_per_iteration"],
                           "phase_name": ["full_run"]*len(full_data_SPR_result["naive_SPR_ll_per_iteration"])
                           }                                            )
        full_iterations_data.to_csv(os.path.join(curr_run_directory, "full_iterations_df.csv"))
        logging.info(f"\n\n Naive SPR  search results: {full_data_SPR_result}")
        return (full_data_SPR_result)
    else:
        all_phases_data = {}
        all_phases_data["lasso_starting_tree_path"] = SPR_chosen_starting_tree_path
        curr_final_tree_object = None
        curr_starting_tree_path = SPR_chosen_starting_tree_path
        all_phases_ll = []
        all_phases_true_ll = []
        all_phases_tag = []
        for i,lasso_configuration in enumerate(lasso_configurations):
            sub_curr_run_directory = os.path.join(curr_run_directory, f"{i}_phase_use_sampled_MSA")
            if not os.path.exists(sub_curr_run_directory):
                os.mkdir(sub_curr_run_directory)
            curr_phase_weights = lasso_configuration["weights_file_path"]
            curr_phase_msa = lasso_configuration["sampled_alignment_path"]
            curr_phase_lasso_intercept = lasso_configuration["lasso_intercept"]
            top_ind = lasso_configuration["top_ind"]
            phase_threshold = lasso_configuration["lasso_threshold"]
            curr_phase_param_dict = SPR_search(
                MSA_path=curr_phase_msa,
                run_unique_name=run_unique_name,
                curr_msa_stats=curr_msa_stats,
                starting_tree_path=curr_starting_tree_path,
                starting_tree_object=curr_final_tree_object,
                curr_run_directory=sub_curr_run_directory,
                top_x_true_trees=top_ind,
                weights_file_path=curr_phase_weights,lasso_intercept = curr_phase_lasso_intercept,n_cpus= curr_msa_stats["n_cpus_Lasso"])
            curr_phase_results_print = {k: curr_phase_param_dict [k] for k in curr_phase_param_dict.keys() if
                                     k not in ["ll_comparison_df", "actual_search_training_df"]
                                     }
            logging.info(f"\n\n{i}'th phase search results: {curr_phase_results_print}")

            curr_phase_param_dict["ll_comparison_df"].to_csv(
                os.path.join(curr_run_directory, "ll_comparison_df_first_phase.csv"))
            actual_search_training_path = os.path.join(curr_msa_stats["curr_msa_version_folder"],
                                                       "actual_search_training_df.csv")
            curr_phase_param_dict["actual_search_training_df"].to_csv(actual_search_training_path)
            prediction_rho_pearson, prediction_pval_pearson, prediction_rho_spearman, prediction_pval_spearman, mse, mistake_cnt = analyze_ll_comparison_df(
                curr_phase_param_dict["ll_comparison_df"])
            ### Continue sampling with full data
            # Use previous tree as a starting tree
            curr_phase_data = {
                f"{i}_phase_lasso_SPR_ll": curr_phase_param_dict["search_best_true_ll"],
                f"{i}_phase_top_ind_tested": top_ind,
                f"{i}_phase_threshold":phase_threshold,
                f"{i}_phase_lasso_SPR_tree_newick": curr_phase_param_dict["search_best_topology_newick"],
                f"{i}_phase_lasso_SPR_spr_moves": curr_phase_param_dict["search_spr_moves"],
                f"{i}_phase_R^2_pearson_during_tree_search": prediction_rho_pearson ** 2,
                f"{i}_phase_R^2_pearson_during_tree_search_pval": prediction_pval_pearson,
                f"{i}_phase_spearmanr_during_tree_search": prediction_rho_spearman,
                f"{i}_phase_spearmanr_during_tree_search_pval": prediction_pval_spearman,
                f"{i}_phase_mse_during_tree_search": mse,
                f"{i}_phase_weights_file_path":curr_phase_weights,
                f"{i}_phase_msa": curr_phase_msa,
                f"{i}_phase_lasso_SPR_starting_tree_path": SPR_chosen_starting_tree_path,
                f"{i}_phase_actual_search_training_path": actual_search_training_path,
                f"{i}_phase_ll_per_iteration": curr_phase_param_dict["ll_per_iteration_list"],
                f"{i}_phase_TRUE_ll_per_iteration": curr_phase_param_dict["TRUE_ll_per_iteration_list"],
                f"{i}_phase_running_time": curr_phase_param_dict["total_search_running_time"],
                f"{i}_phase_no_brlen_running_time": curr_phase_param_dict["total_no_brlen_time"],
                f"{i}_phase_brlen_running_time": curr_phase_param_dict["total_brlen_time"],
                f"{i}_phase_re_opt_running_time": curr_phase_param_dict["total_reoptimization_time"],
                f"{i}_phase_total_spr_neighbours_evaluated": curr_phase_param_dict["total_spr_neighbours_evaluated"],
                f"{i}_phase_total_spr_neighbours_evaluated": curr_phase_param_dict["total_spr_neighbours_evaluated"]

            }
            curr_final_tree_object =  curr_phase_param_dict["search_best_tree_object"]
            curr_starting_tree_path = curr_phase_param_dict["search_best_tree_path"]
            all_phases_ll = all_phases_ll + curr_phase_param_dict["ll_per_iteration_list"]
            all_phases_true_ll = all_phases_true_ll + curr_phase_param_dict["TRUE_ll_per_iteration_list"]
            all_phases_tag =all_phases_tag+ len(curr_phase_param_dict["ll_per_iteration_list"])* [f"{i}_phase"]
            all_phases_data.update(curr_phase_data)

        all_phases_data["best_lasso_tree_newick"] = curr_final_tree_object.write(format=1)
        sub_curr_run_directory = os.path.join(curr_run_directory, f"final_phase_use_sampled_MSA")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
        final_phase_param_dict = SPR_search(
            MSA_path=curr_msa_stats["local_alignment_path"], run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=curr_starting_tree_path,
            starting_tree_object=curr_final_tree_object,
            curr_run_directory=sub_curr_run_directory,
            weights_file_path=False,
            lasso_intercept= -1,
            top_x_true_trees=-1,
            n_cpus= curr_msa_stats["n_cpus_full"]

        )
        final_phase_param_dict["ll_comparison_df"].to_csv(
            os.path.join(curr_run_directory, "ll_comparison_df_final_phase.csv"))

        final_phase_data = {"lasso_SPR_final_phase_ll": final_phase_param_dict["search_best_ll"],
                            "lasso_SPR_final_phase_tree_newick": final_phase_param_dict[
                                "search_best_topology_newick"],
                            "lasso_SPR_final_phase_spr_moves": final_phase_param_dict["search_spr_moves"],
                            "final_phase_ll_per_iteration": final_phase_param_dict["ll_per_iteration_list"],
                            "TRUE_final_phase_ll_per_iteration": final_phase_param_dict[
                                "TRUE_ll_per_iteration_list"],
                            "final_phase_running_time": final_phase_param_dict["total_search_running_time"],
                            "final_phase_no_brlen_running_time": final_phase_param_dict["total_no_brlen_time"],
                            "final_phase_brlen_running_time": final_phase_param_dict["total_brlen_time"],
                            "final_phase_total_spr_neighbours_evaluated": final_phase_param_dict[
                                "total_spr_neighbours_evaluated"]
                            }
        all_phases_ll = all_phases_ll + final_phase_param_dict["ll_per_iteration_list"]
        all_phases_true_ll = all_phases_true_ll + final_phase_param_dict["TRUE_ll_per_iteration_list"]
        all_phases_tag = all_phases_tag + len(final_phase_param_dict["ll_per_iteration_list"]) * ["final_phase"]
        all_phases_data.update(final_phase_data)
        final_optimized_print = {k: final_phase_param_dict[k] for k in final_phase_param_dict.keys()
                                 if
                                 k not in ["ll_comparison_df", "actual_search_training_df"]
                                 }
        logging.info(f"\n\nFinal phase search results: {final_optimized_print}\n")

        all_phases_df = pd.DataFrame({"ll": all_phases_ll, "true_ll": all_phases_true_ll, "phase_name": all_phases_tag})
        all_phases_df.to_csv(os.path.join(curr_run_directory, "all_phases_iterations.csv") )

        return all_phases_data
