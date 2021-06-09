
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


def SPR_iteration(iteration_number,MSA_path, curr_msa_stats, starting_tree_object,
                  curr_run_directory,
                  use_weights):
    add_internal_names(starting_tree_object)
    starting_tree_object.get_tree_root().name="ROOT"
    logging.debug(str(starting_tree_object.write(format=1)) + "\n")
    starting_tree_spr_neighbours =get_possible_spr_moves(get_list_of_edges(starting_tree_object))
    regrafted_trees = [generate_neighbour(starting_tree_object, spr_neighbour) for spr_neighbour in starting_tree_spr_neighbours]
    regrafted_trees_newick = "\n".join([regrafted_tree.write(format=1) for regrafted_tree in regrafted_trees])
    trees_path = os.path.join(curr_run_directory,"iteration_spr_trees")
    with open(trees_path,'w') as TREES:
        TREES.write(regrafted_trees_newick)
    trees_ll, trees_optimized_objects, time_rgft_eval = raxml_optimize_trees_for_given_msa(MSA_path, "rgrft_ll_eval", trees_path,
                                                                           curr_msa_stats, curr_run_directory,
                                                                           weights=curr_msa_stats["weights_file_path"] if use_weights else None)
    if use_weights:
        trees_ll = [(ll/INTEGER_CONST)+curr_msa_stats["lasso_intercept"] for ll in trees_ll]
        trees_true_ll,trees_true_optimized_objects,time_rgft_eval_true = raxml_optimize_trees_for_given_msa(curr_msa_stats["local_alignment_path"],
                                                          "rgrft_ll_eval_on_full_MSA", trees_path,
                                                                                        curr_msa_stats, curr_run_directory,
                                                                                        weights=None)
    else:
        trees_true_ll = trees_ll
    best_ll = max(trees_ll)
    logging.info("Out of {} ll values, the best one found is {}".format(len(trees_ll), best_ll))
    best_ll_index = trees_ll.index(best_ll)
    best_tree_object= trees_optimized_objects[best_ll_index]
    best_true_ll = trees_true_ll[best_ll_index]
    ll_comparison_df = pd.DataFrame(
        {'full msa ll': trees_true_ll, 'sampled msa ll': trees_ll, 'iteration number': iteration_number}
    )
    return best_tree_object, best_ll, best_true_ll,  ll_comparison_df



def get_true_and_local_starting_tree_ll(MSA_path,run_unique_name,starting_tree_path,curr_msa_stats,curr_run_directory,use_weights):
    search_starting_tree_ll,tree_objects, elapsed_running_time_starting_eval = raxml_optimize_trees_for_given_msa(MSA_path,
                                                                       "starting_tree_ll_eval_" + run_unique_name,
                                                                 starting_tree_path,
                                                                 curr_msa_stats,
                                                                 curr_run_directory=curr_run_directory,
                                                                 weights=curr_msa_stats[
                                                                           "weights_file_path"] if use_weights else None)
    if use_weights:
        search_starting_tree_ll = curr_msa_stats["lasso_intercept"]+(search_starting_tree_ll/INTEGER_CONST)
        search_true_starting_tree_ll,tree_objects, elapsed_running_time_starting_eval_true = raxml_optimize_trees_for_given_msa(
            curr_msa_stats["local_alignment_path"], "starting_tree_ll_eval_on_full_" + run_unique_name,
            starting_tree_path,
            curr_msa_stats, curr_run_directory=curr_run_directory, weights=None)
    else:
        search_true_starting_tree_ll = search_starting_tree_ll

    return  search_starting_tree_ll,search_true_starting_tree_ll


def SPR_search(MSA_path, run_unique_name, curr_msa_stats, starting_tree_path, starting_tree_object,
               curr_run_directory,
               use_weights,starting_tree_ll = None ,true_starting_tree_ll = None):
    ll_comparison_df = pd.DataFrame()
    true_vs_sampled_ll_per_iteration_list = []
    spr_iterations_performed_so_far = 0
    if not starting_tree_ll:
        search_starting_tree_ll, search_true_starting_tree_ll = get_true_and_local_starting_tree_ll(MSA_path,run_unique_name,starting_tree_path,curr_msa_stats,curr_run_directory,use_weights)
    else:
        search_starting_tree_ll, search_true_starting_tree_ll = starting_tree_ll, true_starting_tree_ll
    logging.info("Search starting tree ll = {} Search true starting tree ll = {}".format(search_starting_tree_ll, search_true_starting_tree_ll))
    curr_best_tree_ll,curr_best_tree_true_ll = search_starting_tree_ll,search_true_starting_tree_ll
    if starting_tree_object is None:
        curr_best_tree_object = generate_tree_object_from_newick(starting_tree_path)
    else:
        curr_best_tree_object = starting_tree_object
    while True:
        curr_iter_run_directory = os.path.join(curr_run_directory , "iter_" + str(spr_iterations_performed_so_far))
        create_or_clean_dir(curr_iter_run_directory)
        logging.info("iteration number: " + str(spr_iterations_performed_so_far))
        curr_best_neighbour_object, curr_best_neighbour_ll, curr_best_neighbour_true_ll, curr_ll_comparison_df = SPR_iteration(
            spr_iterations_performed_so_far, MSA_path, curr_msa_stats, curr_best_tree_object,
            curr_iter_run_directory,
            use_weights
        )
        logging.info("Our current best tree ll is {} (its true ll is {}), best neighbour ll is {}".format(curr_best_tree_ll, curr_best_tree_true_ll,curr_best_neighbour_ll))
        ll_comparison_df = pd.concat([ll_comparison_df, curr_ll_comparison_df])
        if curr_best_neighbour_ll -  curr_best_tree_ll <= EPSILON:
            logging.info("Difference between best spr neighbour and current tree <= {}, stopping SPR search\n".format(EPSILON))
            break
        logging.info("Updating best neighbour to be our current best tree! ")
        ### Updating current iteration results and preparing for next iteration:
        spr_iterations_performed_so_far = spr_iterations_performed_so_far + 1
        curr_best_tree_object = curr_best_neighbour_object
        curr_best_tree_ll = curr_best_neighbour_ll
        curr_best_tree_true_ll = curr_best_neighbour_true_ll
    return {
        "search_best_ll": curr_best_tree_ll,
        "search_best_true_ll": curr_best_tree_true_ll,
        "search_best_topology_newick": curr_best_tree_object.write(format=1),
        "ll_comparison_df": ll_comparison_df,
        "true_vs_sampled_ll_per_iteration_list": true_vs_sampled_ll_per_iteration_list,
        "search_best_tree_object": curr_best_tree_object,
        "search_spr_moves": spr_iterations_performed_so_far,
    }


def analyze_ll_comparison_df(ll_comparison_df):
    mistake_cnt = 0
    for iteration in ll_comparison_df['iteration number'].unique():
        curr_iter_ll_comparison_df = ll_comparison_df[ll_comparison_df['iteration number'] == iteration]
        maxvalueIndexLabel = curr_iter_ll_comparison_df.idxmax()
        if maxvalueIndexLabel['full msa ll'] != maxvalueIndexLabel['sampled msa ll']:
            mistake_cnt +=1
    if len(ll_comparison_df['full msa ll'])>0:
        rho_pearson, pval_pearson = stats.pearsonr(ll_comparison_df['full msa ll'],
                                                   ll_comparison_df['sampled msa ll'])
        rho_spearman, pval_spearman = stats.spearmanr(ll_comparison_df['full msa ll'],
                                                   ll_comparison_df['sampled msa ll'])
        mse = mean_squared_error(ll_comparison_df['full msa ll'], ll_comparison_df['sampled msa ll'])
    else:
        rho_pearson,pval_pearson,mistake_cnt,rho_spearman, pval_spearman,mse=-1,-1,-1,-1

    return rho_pearson, pval_pearson,rho_spearman, pval_spearman,mse,mistake_cnt







def SPR_analysis(current_file_path,SPR_chosen_starting_tree_path, curr_msa_stats, curr_run_directory,
                 full_run=False):
    run_unique_name="spr"
    logging.info("curr run directory=" + curr_run_directory)
    if not os.path.exists(curr_run_directory):
        os.mkdir(curr_run_directory)
    logging.info(
        "Starting tree is stored in: {}".format(SPR_chosen_starting_tree_path))
    if (full_run):
        logging.info("Starting SPR analysis on full data")
        start_time = time.time()
        full_data_param_dict = SPR_search(
            MSA_path=current_file_path,
            run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=SPR_chosen_starting_tree_path,
            starting_tree_object=None,
            curr_run_directory=curr_run_directory, use_weights = False)
        naive_spr_running_time=time.time()-start_time
        full_data_SPR_result = {"naive_SPR_ll": full_data_param_dict["search_best_ll"],
                                "naive_SPR_spr_moves": full_data_param_dict.get("search_spr_moves"),
                                "naive_SPR_tree_newick" : full_data_param_dict["search_best_topology_newick"]
                                 , "naive_SPR_ll_per_iteration": full_data_param_dict["true_vs_sampled_ll_per_iteration_list"],
                                "naive_SPR_running_time" : naive_spr_running_time
                                }
        logging.info("Full MSA SPR result: " + str(full_data_SPR_result))
        return (full_data_SPR_result)
    else:
        logging.info("Starting SPR analysis on sampled data")
        sub_curr_run_directory = os.path.join(curr_run_directory , "_use_sampled_MSA")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
        start_time = time.time()
        first_optimized_param_dict = SPR_search(
            MSA_path=curr_msa_stats["sampled_alignment_path"],
            run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=SPR_chosen_starting_tree_path,
            starting_tree_object=None,
            curr_run_directory=sub_curr_run_directory,
            use_weights = True)
        first_phase_lasso_running_time = time.time()-start_time
        ll_comparison_df = first_optimized_param_dict["ll_comparison_df"]
        ll_comparison_df.to_csv(os.path.join(curr_run_directory , "ll_comparison_df.csv"))
        prediction_rho_pearson, prediction_pval_pearson,prediction_rho_spearman, prediction_pval_spearman,mse,mistake_cnt = analyze_ll_comparison_df(ll_comparison_df)
        ### Continue sampling with full data
        logging.info("Continue SPR analysis using full data")
        # Use previous tree as a starting tree
        sub_curr_run_directory = os.path.join(curr_run_directory ,"_continue_with_full_MSA")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
        start_time = time.time()
        next_optimized_tree_param_dict = SPR_search(
            MSA_path=curr_msa_stats["local_alignment_path"], run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=None,
            starting_tree_object=first_optimized_param_dict["search_best_tree_object"],
            curr_run_directory=sub_curr_run_directory,
            use_weights = False,
            starting_tree_ll=first_optimized_param_dict["search_best_true_ll"],
            true_starting_tree_ll=first_optimized_param_dict["search_best_true_ll"]
        )
        second_phase_lasso_running_time = time.time() - start_time
        data = {
                "lasso_SPR_first_phase_ll": first_optimized_param_dict["search_best_true_ll"],
                "lasso_SPR_first_phase_tree_newick" : first_optimized_param_dict["search_best_topology_newick"],
                "lasso_SPR_first_phase_spr_moves": first_optimized_param_dict["search_spr_moves"],
                "lasso_SPR_second_phase_ll": next_optimized_tree_param_dict["search_best_ll"],
                "lasso_SPR_second_phase_tree_newick": next_optimized_tree_param_dict["search_best_topology_newick"],
                "lasso_SPR_second_phase_spr_moves": next_optimized_tree_param_dict["search_spr_moves"],
                "R^2_pearson_during_tree_search": prediction_rho_pearson**2,
                "R^2_pearson_during_tree_search_pval": prediction_pval_pearson,
                "spearmanr_during_tree_search": prediction_rho_spearman,
                "spearmanr_during_tree_search_pval": prediction_pval_spearman,
            "first_phase_running_time" : first_phase_lasso_running_time,
            "second_phase_lasso_running_time" : second_phase_lasso_running_time,
                "mse_during_tree_search" : mse,
                "mistake_cnt" : mistake_cnt,
                "lasso_SPR_starting_tree_path": SPR_chosen_starting_tree_path,
            "lasso_ll_per_iteration_first_phase": first_optimized_param_dict["true_vs_sampled_ll_per_iteration_list"],
            "lasso_ll_per_iteration_second_phase": next_optimized_tree_param_dict["true_vs_sampled_ll_per_iteration_list"]

                }
        return data


