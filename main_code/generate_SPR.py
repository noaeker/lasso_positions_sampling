from raxml import *
from config import INTEGER_CONST, EPSILON
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


def get_true_ll_values(curr_msa_stats, trees_path, curr_run_directory, opt_brlen):
    trees_true_ll, trees_true_optimized_objects, time_rgft_eval_true = raxml_optimize_trees_for_given_msa(
        curr_msa_stats["local_alignment_path"],
        "rgrft_ll_eval_on_full_MSA", trees_path,
        curr_msa_stats, curr_run_directory,
        weights=None, n_cpus=curr_msa_stats["n_cpus_full"], opt_brlen=opt_brlen)
    return trees_true_ll


def compute_true_ll_of_best_tree_of_spr_iteration(weights_file_path, best_tree_object, curr_run_directory,
                                                  curr_msa_stats,
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


def regression_correct_lasso_ll_values(lasso_intercept, weights_file_path, trees_ll):
    if weights_file_path:
        if not isinstance(trees_ll, list):
            ll_fixed = (trees_ll + lasso_intercept) / INTEGER_CONST
        else:
            ll_fixed = [((ll) + lasso_intercept) / INTEGER_CONST for ll in trees_ll]
        return ll_fixed
    else:
        return trees_ll


def write_tree_objects_to_file(trees_objects, curr_run_directory, top_trees_file_name):
    if not isinstance(trees_objects, list):
        trees_objects = [trees_objects]
    top_ll_tree_objects = np.array(trees_objects)
    top_ll_trees_newick = "\n".join([tree_object.write(format=1) for tree_object in top_ll_tree_objects])
    top_ll_trees_path = os.path.join(curr_run_directory, top_trees_file_name)
    with open(top_ll_trees_path, 'w') as TOP_LL_TREES:
        TOP_LL_TREES.write(top_ll_trees_newick)
    return top_ll_trees_path


def get_non_greedy_optimized_SPR_neighbours(curr_msa_stats, MSA_path, unique_trees_path, weights_file_path,
                                            lasso_intercept, curr_run_directory, n_cpus):
    if curr_msa_stats["optimized_neighbours_per_iter"] > 1:
        logging.debug("Evaluating (no brlen opt) LL of all SPR neighbours")
        trees_ll_no_brlen, trees_optimized_objects_no_brlen, time_rgft_eval_no_brlen = raxml_optimize_trees_for_given_msa(
            MSA_path, "rgrft_ll_eval_no_brlen",
            unique_trees_path,
            curr_msa_stats,
            curr_run_directory, opt_brlen=False,
            weights=weights_file_path, n_cpus=n_cpus)  # local EVAL of neighbours
        ll_all_spr_candidates_corrected = regression_correct_lasso_ll_values(lasso_intercept, weights_file_path,
                                                                             trees_ll_no_brlen)
        indices_of_spr_candidates_for_brlen_opt = (-np.array(trees_ll_no_brlen)).argsort()[
                                                  :curr_msa_stats["optimized_neighbours_per_iter"]]
        tree_objects_of_spr_candidates_for_brlen_opt = np.array(trees_optimized_objects_no_brlen)[
            indices_of_spr_candidates_for_brlen_opt]
        spr_candidates_for_brlen_opt_file = write_tree_objects_to_file(tree_objects_of_spr_candidates_for_brlen_opt,
                                                                       curr_run_directory,
                                                                       "spr_candidates_for_brlen_opt.trees")
        logging.debug("About to optimize LL of most promising {t} topologies".format(
            t=curr_msa_stats["optimized_neighbours_per_iter"]))
    else:  # No need to eval, performing full branch-lengths optimization on all data
        logging.debug("Fully optimizing all SPR neighbours")
        spr_candidates_for_brlen_opt_file = unique_trees_path
        ll_all_spr_candidates_corrected = []
        time_rgft_eval_no_brlen = 0
    ll_spr_candidates_for_brlen, optimized_objects_spr_candidates_for_brlen, time_rgft_eval_true_spr_candidates_for_brlen = raxml_optimize_trees_for_given_msa(
        MSA_path, "rgrft_ll_eval_brlen",
        spr_candidates_for_brlen_opt_file,
        curr_msa_stats,
        curr_run_directory,
        weights=weights_file_path, n_cpus=n_cpus)
    ll_spr_candidates_for_brlen_corrected = regression_correct_lasso_ll_values(lasso_intercept, weights_file_path,
                                                                               ll_spr_candidates_for_brlen)
    res = {"spr_candidates_for_brlen_opt_file": spr_candidates_for_brlen_opt_file,
           "trees_ll_no_brlen": ll_all_spr_candidates_corrected,
           "ll_spr_candidates_for_brlen_corrected": ll_spr_candidates_for_brlen_corrected,
           "iteration_time_brlen": time_rgft_eval_true_spr_candidates_for_brlen,
           "iteration_time_no_brlen": time_rgft_eval_no_brlen,
           "optimized_objects_spr_candidates_for_brlen": optimized_objects_spr_candidates_for_brlen
           }
    return res


def re_optimize_some_SPR_neighbours_no_weights(ll_spr_candidates_for_brlen_corrected, top_x_true_trees,
                                               optimized_objects_spr_candidates_for_brlen, curr_msa_stats,
                                               curr_run_directory):
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
        weights=None, n_cpus=curr_msa_stats["n_cpus_full"])  # optimize without weights
    best_ll = max(top_trees_true_ll) if isinstance(top_trees_true_ll, list) else top_trees_true_ll
    best_ll_index = top_trees_true_ll.index(best_ll) if isinstance(top_trees_true_ll, list) else 0
    best_tree_object = top_trees_true_optimized_objects[best_ll_index] if isinstance(top_trees_true_optimized_objects,
                                                                                     list) else top_trees_true_optimized_objects

    return best_ll, best_ll_index, best_tree_object, time_rgft_eval_true


# curr_phase_weights = lasso_configuration["weights_file_path"]
# curr_phase_msa = lasso_configuration["sampled_alignment_path"]


def get_first_better_neighbour(prev_lasso_corrected_ll, tree_objects, curr_run_directory, fname, MSA_path,
                               curr_msa_stats, weights_file_path, n_cpus, lasso_intercept, opt_brlen, final_phase=False,
                               final_lasso_configuration=None):
    overall_time = 0
    neighbours_ll = []
    spr_objects = []
    if final_phase:
        logging.debug("final phase, evaluating all neighbours on purpose!")
    for i, candidate_tree in enumerate(tree_objects):
        curr_spr_candidate_for_brlen_opt_file = write_tree_objects_to_file(candidate_tree,
                                                                           curr_run_directory,
                                                                           fname)
        curr_candidate_ll_spr, curr_candidate_object, curr_candidate_object_time = raxml_optimize_trees_for_given_msa(
            MSA_path if not final_phase else final_lasso_configuration["sampled_alignment_path"],
            f"rgrft_ll_eval_{fname}",
            curr_spr_candidate_for_brlen_opt_file,
            curr_msa_stats,
            curr_run_directory,
            weights=final_lasso_configuration["weights_file_path"] if final_phase else weights_file_path, n_cpus=n_cpus,
            opt_brlen=opt_brlen)
        curr_candidate_ll_spr_corrected = regression_correct_lasso_ll_values(lasso_intercept,
                                                                             weights_file_path,
                                                                             curr_candidate_ll_spr)
        overall_time += curr_candidate_object_time
        neighbours_ll.append(curr_candidate_ll_spr_corrected)
        spr_objects.append(curr_candidate_object)
        if (curr_candidate_ll_spr_corrected > prev_lasso_corrected_ll + EPSILON) and not final_phase:
            return curr_candidate_ll_spr_corrected, curr_candidate_object, i + 1, overall_time
    return neighbours_ll, spr_objects, i + 1, overall_time


def find_best_SPR_neighbour_greedy(prev_ll, curr_msa_stats, MSA_path, unique_spr_neighbours_path,
                                   weights_file_path, lasso_intercept,
                                   curr_run_directory, n_cpus, final_phase=False, final_lasso_configuration=None):
    best_tree_object_path = os.path.join(curr_run_directory, "best_greedy_tree")
    logging.debug("Evaluating (greedy) (no brlen opt) LL of SPR neighbours")
    radius_spr_neighbours = generate_multiple_tree_object_from_newick(unique_spr_neighbours_path)
    ll_eval_corrected, spr_objects_eval, n_eval, overall_time_eval = get_first_better_neighbour(prev_ll,
                                                                                                radius_spr_neighbours,
                                                                                                curr_run_directory,
                                                                                                "evaled_trees",
                                                                                                MSA_path,
                                                                                                curr_msa_stats,
                                                                                                weights_file_path,
                                                                                                n_cpus,
                                                                                                lasso_intercept,
                                                                                                opt_brlen=False,
                                                                                                final_phase=final_phase,
                                                                                                final_lasso_configuration=final_lasso_configuration)
    if isinstance(ll_eval_corrected, list):  # i.e., if no better tree was found
        logging.debug(f"No better tree was found using Eval on {n_eval} trees")
        if curr_msa_stats["optimized_neighbours_per_iter"] > -1:
            indices_of_spr_candidates_for_brlen_opt = (-np.array(ll_eval_corrected)).argsort()[
                                                      :curr_msa_stats["optimized_neighbours_per_iter"]]
        else:
            logging.debug("Performing brlen opimization on all trees")
            indices_of_spr_candidates_for_brlen_opt = (-np.array(ll_eval_corrected)).argsort()
        tree_objects_of_spr_candidates_for_brlen_opt = np.array(spr_objects_eval)[
            indices_of_spr_candidates_for_brlen_opt]
        logging.debug("About to optimize (greedy) LL of most promising {treshold} topologies".format(
            treshold=curr_msa_stats["optimized_neighbours_per_iter"]))
        ll_opt_corrected, spr_objects_opt, n_opt, overall_time_opt = get_first_better_neighbour(prev_ll,
                                                                                                tree_objects_of_spr_candidates_for_brlen_opt,
                                                                                                curr_run_directory,
                                                                                                "opt_trees",
                                                                                                MSA_path,
                                                                                                curr_msa_stats,
                                                                                                weights_file_path,
                                                                                                n_cpus,
                                                                                                lasso_intercept,
                                                                                                opt_brlen=True)
        if isinstance(ll_opt_corrected,
                      list):  # if no better tree was found after branch-length optimization, return best obtained tree
            logging.debug(f"No better tree was found after optimizing {n_opt} trees.  ")
            best_opt_ll = max(ll_opt_corrected)
            best_opt_object = spr_objects_opt[ll_opt_corrected.index(best_opt_ll)]
        else:
            best_opt_ll = ll_opt_corrected
            best_opt_object = spr_objects_opt
            logging.debug(
                f" A better tree was found after optimizing all best {n_opt} trees")
        with open(best_tree_object_path, 'w') as GREEDY_TREE:
            GREEDY_TREE.write(best_opt_object.write(format=1))
        tree_opt_true_ll = get_true_ll_values(curr_msa_stats, best_tree_object_path, curr_run_directory, opt_brlen=True)
        logging.debug(
            f"#(prev tree given ll (opt\eval) = ={prev_ll}, curr better tree OPTIMIZED ll= {best_opt_ll}),curr better tree TRUE (FULL DATA) OPTIMIZED ll= {tree_opt_true_ll})")
        return {"best_tree_object": best_opt_object, "best_ll": best_opt_ll, "best_true_ll": tree_opt_true_ll,
                "iteration_time": overall_time_eval + overall_time_opt,
                "iteration_time_brlen": overall_time_opt, "iteration_time_no_brlen": overall_time_eval,
                "n_neighbours_eval": n_eval,
                "n_neighbours_opt": n_opt}

    else:
        best_eval_ll = ll_eval_corrected
        best_eval_object = spr_objects_eval
        logging.debug(f"A better tree was found using Eval on {n_eval} objects. (no need for brlen optimization)")
        with open(best_tree_object_path, 'w') as GREEDY_TREE:
            GREEDY_TREE.write(best_eval_object.write(format=1))
        tree_eval_true_ll = get_true_ll_values(curr_msa_stats, best_tree_object_path, curr_run_directory,
                                               opt_brlen=True)
        logging.debug(
            f"#(prev tree given ll (opt\eval) ={prev_ll}, curr better tree EVAL ll= {best_eval_ll}),curr better tree TRUE (FULL DATA) OPTIMIZED ll= {tree_eval_true_ll})")
        return {"best_tree_object": best_eval_object, "best_ll": best_eval_ll,
                "best_true_ll": tree_eval_true_ll,
                "iteration_time_brlen": 0, "iteration_time_no_brlen": overall_time_eval,
                "n_neighbours_eval": n_eval,
                "n_neighbours_opt": 0, "iteration_time": overall_time_eval}


def find_best_SPR_neighbour_non_greedy(curr_msa_stats, MSA_path, unique_spr_neighbours_path,
                                       weights_file_path, lasso_intercept,
                                       curr_run_directory, iteration_number, top_x_true_trees, n_cpus):
    iteration_time = 0
    re_optimization_time = 0
    spr_evaluation_data = get_non_greedy_optimized_SPR_neighbours(curr_msa_stats, MSA_path, unique_spr_neighbours_path,
                                                                  weights_file_path, lasso_intercept,
                                                                  curr_run_directory, n_cpus)

    iteration_time += spr_evaluation_data["iteration_time_brlen"] + spr_evaluation_data["iteration_time_no_brlen"]
    tree_eval_true_ll = [-1] * len(spr_evaluation_data["trees_ll_no_brlen"])
    tree_opt_true_ll = [-1] * len(spr_evaluation_data["ll_spr_candidates_for_brlen_corrected"])
    if weights_file_path:
        logging.debug(f"SPR iteration {iteration_number} : testing {top_x_true_trees} best Lasso tree objects")
        best_ll, best_ll_index, best_tree_object, re_optimization_time = re_optimize_some_SPR_neighbours_no_weights(
            spr_evaluation_data["ll_spr_candidates_for_brlen_corrected"], top_x_true_trees,
            spr_evaluation_data["optimized_objects_spr_candidates_for_brlen"], curr_msa_stats, curr_run_directory)
        iteration_time += re_optimization_time
        best_true_ll = best_ll
        if curr_msa_stats["compute_all_true_ll"]:
            tree_eval_true_ll = get_true_ll_values(curr_msa_stats, unique_spr_neighbours_path,
                                                   curr_run_directory, opt_brlen=False)
            tree_opt_true_ll = get_true_ll_values(curr_msa_stats,
                                                  spr_evaluation_data["spr_candidates_for_brlen_opt_file"],
                                                  curr_run_directory, opt_brlen=True)
        n_neighbours_reoptimized = top_x_true_trees

    else:
        best_ll = max(spr_evaluation_data["ll_spr_candidates_for_brlen_corrected"])
        best_ll_index = spr_evaluation_data["ll_spr_candidates_for_brlen_corrected"].index(best_ll)
        best_tree_object = spr_evaluation_data["optimized_objects_spr_candidates_for_brlen"][best_ll_index]
        best_true_ll = best_ll
        n_neighbours_reoptimized= 0
    ll_comparison_df_brlen_eval = pd.DataFrame(
        {'full msa ll': tree_eval_true_ll, 'sampled msa ll': spr_evaluation_data["trees_ll_no_brlen"],
         'iteration number': iteration_number}
    )
    ll_comparison_df_brlen_opt = pd.DataFrame(
        {'full msa ll': tree_opt_true_ll,
         'sampled msa ll': spr_evaluation_data["ll_spr_candidates_for_brlen_corrected"],
         'iteration number': iteration_number}
    )
    logging.debug(f"iteration {iteration_number}  best tree ll = {best_ll}")

    results_dict = {"best_tree_object": best_tree_object, "best_ll": best_ll, "best_true_ll": best_true_ll,
                    "ll_comparison_df_brlen_eval": ll_comparison_df_brlen_eval,
                    "ll_comparison_df_brlen_opt": ll_comparison_df_brlen_opt
        , "iteration_time": iteration_time, "iteration_time_brlen": spr_evaluation_data["iteration_time_brlen"],
                    "iteration_time_no_brlen": spr_evaluation_data["iteration_time_no_brlen"],
                    "re_optimization_time": re_optimization_time,
                    "n_neighbours_eval": len(spr_evaluation_data["trees_ll_no_brlen"]),
                    "n_neighbours_opt": len(spr_evaluation_data["ll_spr_candidates_for_brlen_corrected"]),
                    "n_neighbours_reopt":  n_neighbours_reoptimized
                    }
    return results_dict


def SPR_iteration(prev_ll, iteration_number, MSA_path, curr_msa_stats, starting_tree_object,
                  curr_run_directory,
                  weights_file_path, lasso_intercept, top_x_true_trees, n_cpus, final_phase=False,
                  final_lasso_configuration=None):
    add_internal_names(starting_tree_object)
    starting_tree_object.get_tree_root().name = "ROOT"
    logging.debug(str(starting_tree_object.write(format=1)) + "\n")
    starting_tree_spr_neighbours = get_possible_spr_moves(starting_tree_object, rearr_dist=curr_msa_stats["rearr_dist"])
    all_radius_spr_neighbours = [generate_neighbour(starting_tree_object, spr_neighbour) for spr_neighbour in
                                 starting_tree_spr_neighbours]
    regrafted_trees_newick = "\n".join([regrafted_tree.write(format=1) for regrafted_tree in all_radius_spr_neighbours])
    trees_eval_path = os.path.join(curr_run_directory, "iteration_spr_trees")
    with open(trees_eval_path, 'w') as EVAL_TREES:
        EVAL_TREES.write(regrafted_trees_newick)
    unique_spr_neighbours_path = filter_unique_topologies(curr_run_directory, trees_eval_path,
                                                          len(all_radius_spr_neighbours))

    if curr_msa_stats["greedy_SPR"]:
        iteration_results = find_best_SPR_neighbour_greedy(prev_ll, curr_msa_stats, MSA_path,
                                                           unique_spr_neighbours_path,
                                                           weights_file_path, lasso_intercept,
                                                           curr_run_directory, n_cpus, final_phase=final_phase,
                                                           final_lasso_configuration=final_lasso_configuration)
    else:
        iteration_results = find_best_SPR_neighbour_non_greedy(curr_msa_stats, MSA_path, unique_spr_neighbours_path,
                                                               weights_file_path, lasso_intercept,
                                                               curr_run_directory, iteration_number, top_x_true_trees,
                                                               n_cpus)
    return iteration_results


def get_true_and_sampled_starting_tree_ll(reduced_MSA_path, run_unique_name, starting_tree_path, curr_msa_stats,
                                          curr_run_directory, weights_file_path, lasso_intercept, n_cpus):
    search_starting_tree_ll, tree_object_sampled, elapsed_running_time_starting_eval = raxml_optimize_trees_for_given_msa(
        reduced_MSA_path,
        "starting_tree_ll_eval_" + run_unique_name,
        starting_tree_path,
        curr_msa_stats,
        curr_run_directory=curr_run_directory,
        weights=weights_file_path, n_cpus=n_cpus)
    tree_object_true = tree_object_sampled
    if weights_file_path:
        search_starting_tree_ll = \
            regression_correct_lasso_ll_values(lasso_intercept, weights_file_path, [search_starting_tree_ll])[0]
        search_true_starting_tree_ll, tree_object_true, elapsed_running_time_starting_eval_true = raxml_optimize_trees_for_given_msa(
            curr_msa_stats["local_alignment_path"], "starting_tree_ll_eval_on_full_" + run_unique_name,
            starting_tree_path,
            curr_msa_stats, curr_run_directory=curr_run_directory, weights=None, n_cpus=curr_msa_stats["n_cpus_Lasso"])
    else:
        search_true_starting_tree_ll = search_starting_tree_ll

    return search_starting_tree_ll, tree_object_sampled, search_true_starting_tree_ll, tree_object_true


def SPR_search(MSA_path, run_unique_name, curr_msa_stats, starting_tree_path,
               curr_run_directory,
               weights_file_path, lasso_intercept, top_x_true_trees, starting_tree_ll=None, n_cpus=1, final_phase=False,
               final_lasso_configuration=None):
    ll_comparison_df_brlen_eval = pd.DataFrame()
    ll_comparison_df_brlen_opt = pd.DataFrame()
    running_times_per_iter = []
    no_brlen_times_per_iter = []
    brlen_time_per_iter = []
    re_optimization_time_per_iter = []
    spr_neighbours_eval_per_iter = []
    spr_neighbours_opt_per_iter = []
    spr_neighbours_reopt_per_iter = []
    actual_search_training_df = pd.DataFrame()
    spr_iterations_performed_so_far = 0
    search_starting_tree_ll, starting_tree_object_sampled_brlen, search_true_starting_tree_ll, starting_tree_object_true_brlen = get_true_and_sampled_starting_tree_ll(
        MSA_path,
        run_unique_name,
        starting_tree_path,
        curr_msa_stats,
        curr_run_directory,
        weights_file_path,
        lasso_intercept,
        n_cpus=n_cpus)
    LL_per_iteration_list = [search_starting_tree_ll]
    TRUE_LL_per_iteration_list = [search_true_starting_tree_ll]

    logging.info("Search starting tree ll = {} Search true starting tree ll = {}. ".format(search_starting_tree_ll,
                                                                                           search_true_starting_tree_ll))
    if not curr_msa_stats["greedy_SPR"]:
        logging.info("Not greedy search: using true LL values")
        search_starting_tree_ll = search_true_starting_tree_ll
    else:
        logging.info("Greedy search: using local LL values")

    # initialize previous and current trees
    prev_best_tree_ll, prev_best_tree_true_ll = search_starting_tree_ll, search_true_starting_tree_ll
    prev_best_tree_object = starting_tree_object_true_brlen
    while True:
        curr_iter_run_directory = os.path.join(curr_run_directory, "iter_" + str(spr_iterations_performed_so_far))
        create_or_clean_dir(curr_iter_run_directory)
        logging.debug("iteration number: " + str(spr_iterations_performed_so_far))
        new_iteration_results = SPR_iteration(prev_ll=prev_best_tree_ll,
                                              iteration_number=spr_iterations_performed_so_far, MSA_path=MSA_path,
                                              curr_msa_stats=curr_msa_stats,
                                              starting_tree_object=prev_best_tree_object,
                                              curr_run_directory=curr_iter_run_directory,
                                              weights_file_path=weights_file_path, lasso_intercept=lasso_intercept,
                                              top_x_true_trees=top_x_true_trees, n_cpus=n_cpus, final_phase=final_phase,
                                              final_lasso_configuration=final_lasso_configuration
                                              )
        logging.info(
            "Our prev best tree ll is {prev_ll}: (prev true ll = {prev_true}) ; our current best neighbour ll is {best_curr} :(true ll = {best_curr_true})".format(
                prev_ll=prev_best_tree_ll,
                prev_true=prev_best_tree_true_ll, best_curr=new_iteration_results["best_ll"],
                best_curr_true=new_iteration_results["best_true_ll"]

            ))
        ll_comparison_df_brlen_eval = ll_comparison_df_brlen_eval.append(
            new_iteration_results.get("ll_comparison_df_brlen_eval", pd.DataFrame()))
        ll_comparison_df_brlen_opt = ll_comparison_df_brlen_opt.append(
            new_iteration_results.get("ll_comparison_df_brlen_opt", pd.DataFrame()))
        running_times_per_iter.append(new_iteration_results["iteration_time"])
        brlen_time_per_iter.append(new_iteration_results["iteration_time_brlen"])
        no_brlen_times_per_iter.append(new_iteration_results["iteration_time_no_brlen"])
        re_optimization_time_per_iter.append(new_iteration_results.get("re_optimization_time", 0))
        LL_per_iteration_list += [new_iteration_results["best_ll"]]
        TRUE_LL_per_iteration_list += [new_iteration_results["best_true_ll"]]
        spr_neighbours_eval_per_iter.append(new_iteration_results["n_neighbours_eval"])
        spr_neighbours_opt_per_iter.append(new_iteration_results["n_neighbours_opt"])
        spr_neighbours_reopt_per_iter.append(new_iteration_results.get("n_neighbours_reopt", 0))
        spr_iterations_performed_so_far = spr_iterations_performed_so_far + 1

        if new_iteration_results["best_ll"] - prev_best_tree_ll <= EPSILON:  # break if condition holds
            logging.info(
                "Difference between current best spr neighbour and prev tree <= {}, stopping SPR search\n".format(
                    EPSILON))
            break
        ### Updating current iteration results and preparing for next iteration:
        logging.debug("Updating prev neighbour to be our current best tree! ")
        prev_best_tree_ll = new_iteration_results["best_ll"]
        prev_best_tree_true_ll = new_iteration_results["best_true_ll"]
        prev_best_tree_object = new_iteration_results["best_tree_object"]

    search_best_tree_path = os.path.join(curr_run_directory, "search_best_tree_path")
    with open(search_best_tree_path, 'w') as BEST_TREE:
        BEST_TREE.write(prev_best_tree_object.write(format=1))
    search_results = {
        "search_best_ll": prev_best_tree_ll,
        "search_starting_tree_ll": search_starting_tree_ll,
        "search_best_brlen_optimized_true_ll": prev_best_tree_true_ll,
        "search_best_topology_newick": prev_best_tree_object.write(format=1),
        "search_starting_tree_newick": starting_tree_object_true_brlen.write(format=1),
        "ll_comparison_opt": ll_comparison_df_brlen_opt,
        "ll_comparison_eval": ll_comparison_df_brlen_eval,
        "actual_search_training_df": actual_search_training_df,
        "ll_per_iteration_list": LL_per_iteration_list,
        "TRUE_ll_per_iteration_list": TRUE_LL_per_iteration_list,
        "search_best_tree_object": prev_best_tree_object,
        "search_best_tree_path": search_best_tree_path,
        "search_spr_moves": spr_iterations_performed_so_far,
        "search_running_time_per_iter": running_times_per_iter,
        "total_search_running_time": sum(running_times_per_iter),
        "brlen_time_per_iter": brlen_time_per_iter,
        "total_brlen_time": sum(brlen_time_per_iter),
        "no_brlen_time_per_iter": no_brlen_times_per_iter,
        "total_no_brlen_time": sum(no_brlen_times_per_iter),
        "reoptimization_time_per_iter": re_optimization_time_per_iter,
        "total_reoptimization_time": sum(re_optimization_time_per_iter),
        "spr_eval_per_iter": spr_neighbours_eval_per_iter,
        "total_spr_neighbours_eval": sum(spr_neighbours_eval_per_iter),
        "spr_opt_per_iter": spr_neighbours_opt_per_iter,
        "total_spr_neighbours_opt": sum(spr_neighbours_opt_per_iter),
        "spr_reopt_per_iter": spr_neighbours_reopt_per_iter,
        "total_spr_neighbours_reopt": sum(spr_neighbours_reopt_per_iter)
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
        rho_pearson, pval_pearson, mistake_cnt, rho_spearman, pval_spearman, mse = -1, -1, -1, -1, -1, -1

    return rho_pearson, pval_pearson, rho_spearman, pval_spearman, mse, mistake_cnt


def generate_search_correlations_data(param_dict, name,curr_run_directory):
    prediction_rho_pearson_opt, prediction_pval_pearson_opt, prediction_rho_spearman_opt, prediction_pval_spearman_opt, mse_opt, mistake_cnt_opt = analyze_ll_comparison_df(
            param_dict["ll_comparison_opt"])
    prediction_rho_pearson_eval, prediction_pval_pearson_eval, prediction_rho_spearman_eval, prediction_pval_spearman_eval, mse_eval, mistake_cnt_eval = analyze_ll_comparison_df(
            param_dict["ll_comparison_eval"])

    results_dict = {f"{name}_R2_pearson_during_tree_search_eval": prediction_rho_pearson_eval ** 2,
    f"{name}_spearmanr_during_tree_search_eval": prediction_rho_spearman_eval,
    f"{name}_mse_during_tree_search_eval": mse_eval,
                    f"{name}_R2_pearson_during_tree_search_opt": prediction_rho_pearson_opt ** 2,
                    f"{name}_spearmanr_during_tree_search_opt": prediction_rho_spearman_opt,
                    f"{name}_mse_during_tree_search_opt": mse_opt}

    param_dict.get("ll_comparison_opt").to_csv(
        os.path.join(curr_run_directory, f"{name}_ll_comparison_opt.tsv"), index=False, sep='\t')

    param_dict.get("ll_comparison_eval").to_csv(
        os.path.join(curr_run_directory, f"{name}_phase_ll_comparison_eval.tsv"), index=False, sep='\t')
    return results_dict




def generate_per_iter_df(data_param_dict, name):
    iterations_data = pd.DataFrame({"ll": data_param_dict["ll_per_iteration_list"],
                                         "true_ll": data_param_dict["TRUE_ll_per_iteration_list"],
                                         "brlen_times": [0]+data_param_dict["brlen_time_per_iter"],
                                         "no_brlen_times": [0]+data_param_dict["no_brlen_time_per_iter"],
                                        "reopt_times": [0]+ data_param_dict["reoptimization_time_per_iter"],
                                         "total_times": [0]+data_param_dict["search_running_time_per_iter"],
                                         "n_spr_opt": [0]+data_param_dict["spr_opt_per_iter"],
                                         "n_spr_eval": [0]+data_param_dict["spr_eval_per_iter"],
                                         "n_spr_reopt": [0]+data_param_dict["spr_reopt_per_iter"],
                                         "phase_name": [name]*len(data_param_dict["ll_per_iteration_list"])})
    return iterations_data

# 


def generate_search_param_dict(param_dict, name):
    search_final_dict = {
        f"{name}_SPR_ll": param_dict["search_best_brlen_optimized_true_ll"],
        f"{name}_starting_tree_SPR_ll":param_dict["search_starting_tree_ll"],
        f"{name}_best_tree_newick": param_dict["search_best_topology_newick"],
        f"{name}_spr_moves": param_dict["search_spr_moves"],
        f"{name}_running_time": param_dict["total_search_running_time"],
        f"{name}_no_brlen_running_time": param_dict["total_no_brlen_time"],
        f"{name}_brlen_running_time": param_dict["total_brlen_time"],
        f"{name}_re_opt_running_time": param_dict["total_reoptimization_time"],
        f"{name}_total_spr_neighbours_evaluated": param_dict["total_spr_neighbours_eval"],
        f"{name}_total_spr_neighbours_optimized": param_dict["total_spr_neighbours_opt"],
        f"{name}_total_spr_neighbours_reoptimized": param_dict["total_spr_neighbours_reopt"]
    }
    return search_final_dict


def SPR_analysis(lasso_configurations, SPR_chosen_starting_tree_path, curr_msa_stats, curr_run_directory,
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
            curr_run_directory=curr_run_directory, weights_file_path=None, lasso_intercept=-1, top_x_true_trees=1,
            n_cpus=curr_msa_stats["n_cpus_full"])
        full_data_SPR_result = generate_search_param_dict(full_data_param_dict, "naive_SPR")
        full_iterations_data = generate_per_iter_df( full_data_param_dict, "full_data")
        full_iterations_data.to_csv(os.path.join(curr_run_directory, "full_iterations_df.tsv"), index=False, sep='\t')
        naive_search_resluts_print = {k: full_data_SPR_result[k] for k in full_data_SPR_result.keys() if
                                      all(x not in k for x in ["path", "newick"])
                                      }
        logging.info(f"\n\n Naive SPR  search results: {naive_search_resluts_print}")
        return (full_data_SPR_result)
    else:
        all_phases_data = {}
        all_phases_iterations_data = pd.DataFrame()
        all_phases_data["lasso_starting_tree_path"] = SPR_chosen_starting_tree_path
        curr_final_tree_object = None
        curr_starting_tree_path = SPR_chosen_starting_tree_path
        for i, lasso_configuration in enumerate(lasso_configurations):
            phase_name = f"phase_{i}"
            sub_curr_run_directory = os.path.join(curr_run_directory, f"{phase_name}_use_sampled_MSA")
            if not os.path.exists(sub_curr_run_directory):
                os.mkdir(sub_curr_run_directory)
            curr_phase_weights = lasso_configuration["weights_file_path"]
            curr_phase_msa_path = lasso_configuration["sampled_alignment_path"]
            curr_phase_lasso_intercept = lasso_configuration["lasso_intercept"]
            top_ind = lasso_configuration["top_ind"]
            phase_threshold = lasso_configuration["lasso_threshold"]
            curr_phase_param_dict = SPR_search(
                MSA_path=curr_phase_msa_path,
                run_unique_name=run_unique_name,
                curr_msa_stats=curr_msa_stats,
                starting_tree_path=curr_starting_tree_path,
                curr_run_directory=sub_curr_run_directory,
                top_x_true_trees=top_ind,
                weights_file_path=curr_phase_weights, lasso_intercept=curr_phase_lasso_intercept,
                n_cpus=curr_msa_stats["n_cpus_Lasso"])

            curr_phase_data = generate_search_param_dict(curr_phase_param_dict, name = phase_name)
            curr_phase_data.update({
            f"{phase_name}_top_ind_tested": top_ind,
            f"{phase_name}_threshold": phase_threshold,
            f"{phase_name}_weights_file_path": curr_phase_weights,
            f"{phase_name}_msa_path": curr_phase_msa_path,
            f"{phase_name}_lasso_SPR_starting_tree_path": SPR_chosen_starting_tree_path})

            if not curr_msa_stats["greedy_SPR"] and curr_msa_stats["compute_all_true_ll"]:
                curr_phase_data.update(generate_search_correlations_data(curr_msa_stats,curr_phase_param_dict,phase_name,curr_run_directory = sub_curr_run_directory))

            curr_phase_results_print = {k: curr_phase_data[k] for k in curr_phase_data.keys() if
                                        all(x not in k for x in ["path", "pval", "newick"])
                                        }
            logging.info(f"\n\n{i}'th phase search results: {curr_phase_results_print}")

            curr_final_tree_object = curr_phase_param_dict["search_best_tree_object"]
            curr_starting_tree_path = curr_phase_param_dict["search_best_tree_path"]
            curr_phase_iterations_data = generate_per_iter_df(curr_phase_param_dict, f"{i}_phase")
            all_phases_iterations_data = all_phases_iterations_data.append(curr_phase_iterations_data)
            all_phases_data.update(curr_phase_data)

        all_phases_data["best_lasso_tree_newick"] = curr_final_tree_object.write(format=1)
        sub_curr_run_directory = os.path.join(curr_run_directory, f"final_phase_use_sampled_MSA")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
        if curr_msa_stats["use_modified_final_search"]:
            final_phase_param_dict = SPR_search(
                MSA_path=curr_msa_stats["local_alignment_path"], run_unique_name=run_unique_name,
                curr_msa_stats=curr_msa_stats,
                starting_tree_path=curr_starting_tree_path,
                curr_run_directory=sub_curr_run_directory,
                weights_file_path=False,
                lasso_intercept=-1,
                top_x_true_trees=-1,
                n_cpus=curr_msa_stats["n_cpus_full"],
                final_phase=True,
                final_lasso_configuration=lasso_configurations[-1]

        )
        else:
            final_phase_param_dict = SPR_search(
                MSA_path=curr_msa_stats["local_alignment_path"], run_unique_name=run_unique_name,
                curr_msa_stats=curr_msa_stats,
                starting_tree_path=curr_starting_tree_path,
                curr_run_directory=sub_curr_run_directory,
                weights_file_path=False,
                lasso_intercept=-1,
                top_x_true_trees=-1,
                n_cpus=curr_msa_stats["n_cpus_full"]

        )
        final_phase_data = generate_search_param_dict(final_phase_param_dict, name="final_phase")
        all_phases_data.update(final_phase_data)
        final_optimized_print = {k: final_phase_data[k] for k in final_phase_data.keys()
                                 if
                                 all(x not in k for x in ["path", "pval", "newick"])
                                 }
        logging.info(f"\n\nFinal phase search results: {final_optimized_print}\n")

        final_phase_iterations_data = generate_per_iter_df(final_phase_param_dict, "final_phase")
        all_phases_iterations_data = all_phases_iterations_data.append(final_phase_iterations_data)
        all_phases_iterations_data.to_csv(os.path.join(curr_run_directory, "all_phases_iterations.tsv"), index=False, sep='\t')

        return all_phases_data
