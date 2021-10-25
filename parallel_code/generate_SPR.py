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


def get_true_ll_values_and_sitelh(use_weights, curr_msa_stats, trees_path, curr_run_directory, trees_ll):
    if use_weights:
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
                weights=None)
        return trees_true_ll, pd.DataFrame()
    else:
        return trees_ll, pd.DataFrame()


def compute_true_ll_of_best_tree_of_spr_iteration(use_weights, best_tree_object, curr_run_directory, curr_msa_stats,
                                                  trees_true_ll, best_ll_index, best_ll, top_x_to_test):
    if use_weights and top_x_to_test == 1:
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
                weights=None)
            best_true_ll = best_tree_true_ll
            return best_true_ll
    else:
        return best_ll


def regression_correct_lasso_ll_values(use_weights, curr_msa_stats, trees_ll):
    if use_weights:
        ll_fixed = [((ll) + curr_msa_stats["lasso_intercept"]) / INTEGER_CONST for ll in trees_ll]
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





def get_locally_optimized_SPR_neighbours(curr_msa_stats, MSA_path, unique_trees_path, use_weights, curr_run_directory):
    if curr_msa_stats["optimized_neighbours_per_iter"] > 1:
        logging.info("Evaluating (no brlen opt) LL of all SPR neighbours")
        trees_ll_no_brlen, trees_optimized_objects_no_brlen, time_rgft_eval_no_brlen = raxml_optimize_trees_for_given_msa(
            MSA_path, "rgrft_ll_eval_no_brlen",
            unique_trees_path,
            curr_msa_stats,
            curr_run_directory, opt_brlen=False,
            weights=curr_msa_stats[
                "weights_file_path"] if use_weights else None)
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
        weights=curr_msa_stats[
            "weights_file_path"] if use_weights else None)
    ll_spr_candidates_for_brlen_corrected = regression_correct_lasso_ll_values(use_weights, curr_msa_stats,
                                                                               ll_spr_candidates_for_brlen)
    basic_neighbours_optimization_time = time_rgft_eval_true_spr_candidates_for_brlen+time_rgft_eval_no_brlen
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
        weights=None)  # optimize without weights
    best_ll = max(top_trees_true_ll)
    best_ll_index = top_trees_true_ll.index(best_ll)
    best_tree_object = top_trees_true_optimized_objects[best_ll_index]
    return best_ll, best_ll_index, best_tree_object,time_rgft_eval_true



def SPR_iteration(iteration_number, MSA_path, curr_msa_stats, starting_tree_object,
                  curr_run_directory,
                  use_weights, top_x_true_trees):
    iteration_time = 0
    re_optimization_time=-1
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
    ll_spr_candidates_for_brlen_corrected, optimized_objects_spr_candidates_for_brlen, time_rgft_eval_true_spr_candidates_for_brlen,time_rgft_eval_no_brlen = get_locally_optimized_SPR_neighbours(curr_msa_stats, MSA_path, unique_trees_path, use_weights, curr_run_directory)
    iteration_time+= time_rgft_eval_true_spr_candidates_for_brlen+time_rgft_eval_no_brlen
    if top_x_true_trees > 1 and use_weights:
        logging.debug(f"SPR iteration {iteration_number} : testing {top_x_true_trees} best Lasso tree objects")
        best_ll, best_ll_index, best_tree_object, re_optimization_time = re_optimize_some_SPR_neighbours_no_weights(ll_spr_candidates_for_brlen_corrected,top_x_true_trees,optimized_objects_spr_candidates_for_brlen, curr_msa_stats, curr_run_directory)
        iteration_time += re_optimization_time
    else:
        best_ll = max(ll_spr_candidates_for_brlen_corrected)
        best_ll_index = ll_spr_candidates_for_brlen_corrected.index(best_ll)
        best_tree_object = optimized_objects_spr_candidates_for_brlen[best_ll_index]
    trees_true_ll, true_sitelh_df = get_true_ll_values_and_sitelh(use_weights, curr_msa_stats, trees_path,
                                                                  curr_run_directory, ll_spr_candidates_for_brlen_corrected)
    best_true_ll = compute_true_ll_of_best_tree_of_spr_iteration(use_weights, best_tree_object, curr_run_directory,
                                                                 curr_msa_stats, trees_true_ll, best_ll_index, best_ll,
                                                                 top_x_true_trees)
    ll_comparison_df = pd.DataFrame(
        {'full msa ll': trees_true_ll, 'sampled msa ll': ll_spr_candidates_for_brlen_corrected, 'iteration number': iteration_number}
    )

    return best_tree_object, best_ll, best_true_ll, ll_comparison_df, true_sitelh_df, iteration_time,time_rgft_eval_true_spr_candidates_for_brlen,time_rgft_eval_no_brlen,re_optimization_time



def get_true_and_sampled_starting_tree_ll(MSA_path, run_unique_name, starting_tree_path, curr_msa_stats,
                                          curr_run_directory, use_weights):
    search_starting_tree_ll, tree_objects, elapsed_running_time_starting_eval = raxml_optimize_trees_for_given_msa(
        MSA_path,
        "starting_tree_ll_eval_" + run_unique_name,
        starting_tree_path,
        curr_msa_stats,
        curr_run_directory=curr_run_directory,
        weights=curr_msa_stats[
            "weights_file_path"] if use_weights else None)
    if use_weights:
        search_starting_tree_ll = regression_correct_lasso_ll_values(True, curr_msa_stats, [search_starting_tree_ll])[0]
        search_true_starting_tree_ll, tree_objects, elapsed_running_time_starting_eval_true = raxml_optimize_trees_for_given_msa(
            curr_msa_stats["local_alignment_path"], "starting_tree_ll_eval_on_full_" + run_unique_name,
            starting_tree_path,
            curr_msa_stats, curr_run_directory=curr_run_directory, weights=None)
    else:
        search_true_starting_tree_ll = search_starting_tree_ll

    return search_starting_tree_ll, search_true_starting_tree_ll


def SPR_search(MSA_path, run_unique_name, curr_msa_stats, starting_tree_path, starting_tree_object,
               curr_run_directory,
               use_weights, top_x_true_trees, starting_tree_ll=None, true_starting_tree_ll=None):
    ll_comparison_df = pd.DataFrame()
    running_times_per_iter = []
    no_brlen_times_per_iter = []
    brlen_per_iter = []
    re_optimization_time_per_iter = []
    actual_search_training_df = pd.DataFrame()
    spr_iterations_performed_so_far = 0
    if not starting_tree_ll:
        search_starting_tree_ll, search_true_starting_tree_ll = get_true_and_sampled_starting_tree_ll(MSA_path,
                                                                                                      run_unique_name,
                                                                                                      starting_tree_path,
                                                                                                      curr_msa_stats,
                                                                                                      curr_run_directory,
                                                                                                      use_weights)

    else:
        search_starting_tree_ll, search_true_starting_tree_ll = starting_tree_ll, true_starting_tree_ll
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
        curr_best_neighbour_object, curr_best_neighbour_ll, curr_best_neighbour_true_ll, curr_ll_comparison_df, curr_true_sitelh_df, iteration_time,time_rgft_eval_true_spr_candidates_for_brlen,time_rgft_eval_no_brlen,re_optimization_time = SPR_iteration(
            spr_iterations_performed_so_far, MSA_path, curr_msa_stats, curr_best_tree_object,
            curr_iter_run_directory,
            use_weights,
            top_x_true_trees
        )
        logging.debug(
            "Our current best tree ll is {} (its true ll is {}), best neighbour ll is {}".format(curr_best_tree_ll,
                                                                                                 curr_best_tree_true_ll,
                                                                                                 curr_best_neighbour_ll))
        ll_comparison_df = ll_comparison_df.append(curr_ll_comparison_df)
        running_times_per_iter.append(iteration_time)
        brlen_per_iter.append(time_rgft_eval_true_spr_candidates_for_brlen)
        no_brlen_times_per_iter.append(time_rgft_eval_no_brlen)
        re_optimization_time_per_iter.append(re_optimization_time)

        actual_search_training_df = actual_search_training_df.append(curr_true_sitelh_df)
        if curr_best_neighbour_ll - curr_best_tree_ll <= EPSILON:
            logging.debug(
                "Difference between best spr neighbour and current tree <= {}, stopping SPR search\n".format(EPSILON))
            break
        logging.debug("Updating best neighbour to be our current best tree! ")
        ### Updating current iteration results and preparing for next iteration:
        spr_iterations_performed_so_far = spr_iterations_performed_so_far + 1
        curr_best_tree_object = curr_best_neighbour_object
        curr_best_tree_ll = curr_best_neighbour_ll
        curr_best_tree_true_ll = curr_best_neighbour_true_ll
        LL_per_iteration_list += [curr_best_tree_ll]
        TRUE_LL_per_iteration_list += [curr_best_tree_true_ll]

    search_results = {
        "search_best_ll": curr_best_tree_ll,
        "search_starting_tree_ll": starting_tree_ll,
        "search_best_true_ll": curr_best_tree_true_ll,
        "search_best_topology_newick": curr_best_tree_object.write(format=1),
        "search_starting_tree_newick": starting_tree_object.write(format=1),
        "ll_comparison_df": ll_comparison_df,
        "actual_search_training_df": curr_true_sitelh_df,
        "ll_per_iteration_list": LL_per_iteration_list,
        "TRUE_ll_per_iteration_list": TRUE_LL_per_iteration_list,
        "search_best_tree_object": curr_best_tree_object,
        "search_spr_moves": spr_iterations_performed_so_far,
        "running_time_per_iter": running_times_per_iter,
        "total_search_running_time": sum(running_times_per_iter),
        "total_brlen_time" : sum(brlen_per_iter) ,
        "total_no_brlen_time" : sum(no_brlen_times_per_iter),
        "total_reoptimization_time" : sum (re_optimization_time_per_iter)
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


def SPR_analysis(current_file_path, SPR_chosen_starting_tree_path, curr_msa_stats, curr_run_directory,
                 full_run=False):
    run_unique_name = "spr"
    if not os.path.exists(curr_run_directory):
        os.mkdir(curr_run_directory)
    if (full_run):
        full_data_param_dict = SPR_search(
            MSA_path=current_file_path,
            run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=SPR_chosen_starting_tree_path,
            starting_tree_object=None,
            curr_run_directory=curr_run_directory, use_weights=False, top_x_true_trees=1)
        full_data_SPR_result = {"naive_SPR_ll": full_data_param_dict["search_best_ll"],
                                "naive_SPR_spr_moves": full_data_param_dict["search_spr_moves"],
                                "naive_SPR_tree_newick": full_data_param_dict["search_best_topology_newick"]
            , "naive_SPR_ll_per_iteration": full_data_param_dict["ll_per_iteration_list"],
                                "SPR_search_starting_tree_ll": full_data_param_dict["search_starting_tree_ll"],
                                "SPR_search_starting_tree_newick": full_data_param_dict["search_starting_tree_newick"],
                                "naive_spr_brlen_running_time" : full_data_param_dict["total_brlen_time"],
                                "naive_spr_no_brlen_running_time": full_data_param_dict["total_no_brlen_time"],
                                "naive_SPR_running_time": full_data_param_dict["total_search_running_time"]


                                }
        full_iterations_data = pd.DataFrame({"ll": full_data_SPR_result["naive_SPR_ll_per_iteration"],"true_ll": full_data_SPR_result["naive_SPR_ll_per_iteration"],
                           "phase_name": ["full_run"]*len(full_data_SPR_result["naive_SPR_ll_per_iteration"])
                           }                                            )
        full_iterations_data.to_csv(os.path.join(curr_run_directory, "full_iterations_df.csv"))
        logging.info(f"\n\nFinal phase search results: {full_data_SPR_result}")
        return (full_data_SPR_result)
    else:
        sub_curr_run_directory = os.path.join(curr_run_directory, "_use_sampled_MSA_first_phase")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
        first_phase_param_dict = SPR_search(
            MSA_path=curr_msa_stats["sampled_alignment_path"],
            run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=SPR_chosen_starting_tree_path,
            starting_tree_object=None,
            curr_run_directory=sub_curr_run_directory,
            top_x_true_trees=curr_msa_stats["top_ind_to_test_first_phase"],
            use_weights=True)
        first_optimized_print = {k: first_phase_param_dict[k] for k in first_phase_param_dict.keys() if
                                 k not in ["ll_comparison_df", "actual_search_training_df"]
                                 }
        logging.info(f"\n\nFirst phase search results: {first_optimized_print}")

        first_phase_param_dict["ll_comparison_df"].to_csv(
            os.path.join(curr_run_directory, "ll_comparison_df_first_phase.csv"))
        actual_search_training_path = os.path.join(curr_msa_stats["curr_msa_version_folder"],
                                                   "actual_search_training_df.csv")
        first_phase_param_dict["actual_search_training_df"].to_csv(actual_search_training_path)
        prediction_rho_pearson, prediction_pval_pearson, prediction_rho_spearman, prediction_pval_spearman, mse, mistake_cnt = analyze_ll_comparison_df(
            first_phase_param_dict["ll_comparison_df"])
        ### Continue sampling with full data
        # Use previous tree as a starting tree
        first_phase_data = {
            "lasso_SPR_first_phase_ll": first_phase_param_dict["search_best_true_ll"],
            "lasso_SPR_first_phase_tree_newick": first_phase_param_dict["search_best_topology_newick"],
            "lasso_SPR_first_phase_spr_moves": first_phase_param_dict["search_spr_moves"],
            "R^2_pearson_during_tree_search": prediction_rho_pearson ** 2,
            "R^2_pearson_during_tree_search_pval": prediction_pval_pearson,
            "spearmanr_during_tree_search": prediction_rho_spearman,
            "spearmanr_during_tree_search_pval": prediction_pval_spearman,
            "mse_during_tree_search": mse,
            "lasso_SPR_starting_tree_path": SPR_chosen_starting_tree_path,
            "actual_search_training_path": actual_search_training_path,
            "first_phase_ll_per_iteration": first_phase_param_dict["ll_per_iteration_list"],
            "TRUE_first_phase_ll_per_iteration": first_phase_param_dict["TRUE_ll_per_iteration_list"],
            "first_phase_running_time": first_phase_param_dict["total_search_running_time"],
            "first_phase_no_brlen_running_time": first_phase_param_dict["total_no_brlen_time"],
            "first_phase_brlen_running_time": first_phase_param_dict["total_brlen_time"],
            "first_phase_re_opt_running_time": first_phase_param_dict["total_reoptimization_time"]

        }
        sub_curr_run_directory = os.path.join(curr_run_directory, "_use_sampled_MSA_second_phase")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
        second_phase_param_dict = SPR_search(
            MSA_path=curr_msa_stats["sampled_alignment_path"], run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=None,
            starting_tree_object=first_phase_param_dict["search_best_tree_object"],
            curr_run_directory=sub_curr_run_directory,
            use_weights=True,
            top_x_true_trees=curr_msa_stats["top_ind_to_test_second_phase"],
            starting_tree_ll=first_phase_param_dict["search_best_ll"],
            true_starting_tree_ll=first_phase_param_dict["search_best_true_ll"]

        )
        second_phase_param_dict["ll_comparison_df"].to_csv(
            os.path.join(curr_run_directory, "ll_comparison_df_second_phase.csv"))

        second_phase_data = {"lasso_SPR_second_phase_ll": second_phase_param_dict["search_best_ll"],
                             "lasso_SPR_second_phase_tree_newick": second_phase_param_dict[
                                 "search_best_topology_newick"],
                             "lasso_SPR_second_phase_spr_moves": second_phase_param_dict["search_spr_moves"],
                             "second_phase_ll_per_iteration": second_phase_param_dict["ll_per_iteration_list"],
                             "TRUE_second_phase_ll_per_iteration": second_phase_param_dict["TRUE_ll_per_iteration_list"],
                             "second_phase_running_time": second_phase_param_dict["total_search_running_time"],
                             "second_phase_no_brlen_running_time": second_phase_param_dict["total_no_brlen_time"],
                             "second_phase_brlen_running_time": second_phase_param_dict["total_brlen_time"],
                             "second_phase_re_opt_running_time": second_phase_param_dict["total_reoptimization_time"]
                             }
        second_optimized_print = {k: second_phase_param_dict[k] for k in second_phase_param_dict.keys() if
                                  k not in ["ll_comparison_df", "actual_search_training_df"]
                                  }
        logging.info(f"\n\nSecond phase search results: {second_optimized_print}")

        sub_curr_run_directory = os.path.join(curr_run_directory, "_finalize_with_full_MSA")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)

        final_phase_param_dict = SPR_search(
            MSA_path=curr_msa_stats["local_alignment_path"], run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=None,
            starting_tree_object=second_phase_param_dict["search_best_tree_object"],
            curr_run_directory=sub_curr_run_directory,
            use_weights=False,
            top_x_true_trees=-1,
            starting_tree_ll=second_phase_param_dict["search_best_true_ll"],
            true_starting_tree_ll=second_phase_param_dict["search_best_true_ll"]

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
                            "final_phase_brlen_running_time": final_phase_param_dict["total_brlen_time"]
                            }
        final_optimized_print = {k: final_phase_param_dict[k] for k in final_phase_param_dict.keys()
                                 if
                                 k not in ["ll_comparison_df", "actual_search_training_df"]
                                 }
        logging.info(f"\n\nFinal phase search results: {final_optimized_print}\n")
        first_phase_data.update(second_phase_data)
        first_phase_data.update(final_phase_data)

        all_phases_ll = first_phase_data["first_phase_ll_per_iteration"] + first_phase_data[
            "second_phase_ll_per_iteration"] + first_phase_data["final_phase_ll_per_iteration"]
        all_phases_true_ll = first_phase_data["TRUE_first_phase_ll_per_iteration"] + first_phase_data[
            "TRUE_second_phase_ll_per_iteration"] + first_phase_data["TRUE_final_phase_ll_per_iteration"]
        all_phases = len(first_phase_data["first_phase_ll_per_iteration"]) * ["first_phase"] + len(
            second_phase_data["second_phase_ll_per_iteration"]) * ["second_phase"] + len(
            final_phase_data["final_phase_ll_per_iteration"]) * ["final_phase" ]

        all_phases_df = pd.DataFrame({"ll": all_phases_ll, "true_ll": all_phases_true_ll, "phase_name": all_phases})
        all_phases_df.to_csv(os.path.join(curr_run_directory, "all_phases_iterations.csv") )

        return first_phase_data
