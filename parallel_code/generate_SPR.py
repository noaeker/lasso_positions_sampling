
from raxml import *
import numpy as np
import shutil
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


def SPR_iteration(MSA_path, curr_msa_stats, starting_tree_object, starting_tree_ll, starting_tree_ll_on_full_data,
                  iteration_number, curr_run_directory,
                  best_topology_path, spr_log_file_object, use_weights):
    ll_comparison_df = pd.DataFrame(columns=['full msa ll', 'sampled msa ll', 'iteration number'])
    starting_tree = starting_tree_object
    best_tree_object = starting_tree_object
    logging.debug("\n\niteration starting tree=\n" + starting_tree_object.get_ascii(attributes=['name'],
                                                                                    show_internal=True) + "\n")
    logging.debug(str(starting_tree.write(format=1)) + "\n")
    best_ll = starting_tree_ll
    best_true_ll = starting_tree_ll_on_full_data
    neighbours_tested = 0
    for i, prune_node in enumerate(starting_tree.iter_descendants("levelorder")):
        pruned_node_parent_name, pruned_subtree, remaining_tree = prune_at_internal_node(starting_tree,
                                                                                         prune_node.name)  # subtree1 is the pruned subtree. subtree2 is the remaining subtree
        for j, rgft_node in enumerate(remaining_tree.iter_descendants("levelorder")):
            if prune_node.up.name != rgft_node.up.name:
                if DETAILED_SPR_LOG:
                    spr_log_file_object.write(
                        "prune node parent name= {} rgft node parent name {}".format(prune_node.up.name,
                                                                                     rgft_node.up.name) + "\n")
                logging.debug("i=" + str(i) + " j=" + str(j))
                if DETAILED_SPR_LOG:
                    spr_log_file_object.write("j = {} ; Refraft node name = {}".format(j, rgft_node.name) + "\n")
                rgrft_folder = os.path.join(curr_run_directory , "_i_" + str(i) + "_j_" + str(
                    j) + "_iter_" + str(iteration_number))
                create_or_clean_dir(rgrft_folder)
                rgrft_path = os.path.join(rgrft_folder , rgft_node.name)
                if not os.path.exists(rgrft_path):
                    regrafted_tree = regraft_as_sister_of_given_internal_node(rgft_node.name,
                                                                              pruned_subtree, remaining_tree)
                    if DETAILED_SPR_LOG:
                        print_subtree(pruned_subtree, spr_log_file_object, text="pruned subtree")
                        print_subtree(remaining_tree, spr_log_file_object, text="Remaining tree")
                        print_subtree(regrafted_tree, spr_log_file_object, text="regrafted tree ")
                    regrafted_tree.write(format=1, outfile=rgrft_path)
                    logging.debug("evaluating regrafted ll on given data in : " + MSA_path)
                    ll = raxml_optimize_ll_on_given_tree_and_msa(MSA_path, "rgrft_ll_eval", rgrft_path,
                                                                 curr_msa_stats, rgrft_folder,
                                                                 weights=curr_msa_stats["weights_file_path"] if use_weights else None)
                    logging.debug("evaluating regrafted ll on true data in : " + curr_msa_stats["local_alignment_path"])
                    true_ll = raxml_optimize_ll_on_given_tree_and_msa(curr_msa_stats["local_alignment_path"],
                                                                          "rgrft_ll_eval_on_full_MSA", rgrft_path,
                                                                      curr_msa_stats, rgrft_folder,
                                                                      weights=None)
                    logging.debug("curr regraft ll=" + str(ll) + " and curr true ll= " + str(true_ll))
                    neighbours_tested = neighbours_tested + 1

                    logging.debug(
                        "Current SPR neighboor log likelihood={}, best log likelihood = {}".format(ll, best_ll))
                    ll_comparison_df = ll_comparison_df.append(
                        {'full msa ll': true_ll, 'sampled msa ll': ll, 'iteration number': iteration_number},
                        ignore_index=True
                    )
                if DETAILED_SPR_LOG:
                    write_spr_log_message(spr_log_file_object, rgrft_path, best_ll, ll, best_topology_path)
                if ll > best_ll:  # updating best object
                    logging.debug("updating best ll to be " + str(ll) + " and best true ll to be " + str(true_ll))
                    best_ll = ll
                    best_true_ll = true_ll
                    best_tree_object = regrafted_tree
                    shutil.copy(rgrft_path, best_topology_path)
                if os.path.exists(rgrft_folder) and DELETE_SPR_FILES:
                    logging.debug("Removing rgrft folder:" + rgrft_folder)
                    if DETAILED_SPR_LOG:
                        spr_log_file_object.write("Removing dir:" + rgrft_folder + "\n")
                    shutil.rmtree(rgrft_folder)
    if DETAILED_SPR_LOG:
        print_subtree(best_tree_object, spr_log_file_object, text="********iteration best tree")
    # print_subtree(best_tree_object, None, text="********iteration best tree")
    # logging.info("Number of SPR neighbours tested is: " + str(neighbours_tested))
    # logging.info("iteration best ll =" + str(best_ll) + " and iteration best true ll=" + str(best_true_ll))

    return best_tree_object, best_ll, best_true_ll, neighbours_tested, ll_comparison_df










def SPR_search(MSA_path, run_unique_name, curr_msa_stats, starting_tree_path, starting_tree_object,
               curr_run_directory,
               phase_name,
               use_weights):
    ll_comparison_df = pd.DataFrame(columns=['full msa ll', 'sampled msa ll', 'iteration number'])
    spr_log_file = curr_msa_stats["spr_log_path"]
    with open(spr_log_file, 'a') as spr_log_file_object:
        number_of_SPR_neighbours_per_iteration_list = []
        true_vs_sampled_ll_per_iteration_list = []
        spr_log_file_object.write("\n\n\nStarting new SPR search on :" + run_unique_name + "\n")
        # logging.info("Starting naive SPR search on file: " + MSA_path)
        spr_iterations_performed_so_far = 0
        curr_iter_starting_tree_topology = starting_tree_path
        # logging.info("evaluating starting tree on curr MSA: " + MSA_path)
        naive_SPR_search_starting_tree_ll = raxml_optimize_ll_on_given_tree_and_msa(MSA_path,
                                                                                        "starting_tree_ll_eval_" + run_unique_name,
                                                                                    starting_tree_path,
                                                                                    curr_msa_stats,
                                                                                    curr_run_directory=curr_run_directory,
                                                                                    weights=curr_msa_stats["weights_file_path"] if use_weights else None)
        # logging.info("evaluating starting tree on full MSA: " + curr_msa_stats["local_alignment_path"])
        naive_SPR_search_starting_tree_true_ll = raxml_optimize_ll_on_given_tree_and_msa(
            curr_msa_stats["local_alignment_path"], "starting_tree_ll_eval_on_full_" + run_unique_name,
            starting_tree_path,
            curr_msa_stats, curr_run_directory=curr_run_directory, weights=None)
        curr_iter_starting_tree_ll = naive_SPR_search_starting_tree_ll
        curr_iter_starting_tree_ll_on_full_data = naive_SPR_search_starting_tree_true_ll
        if starting_tree_object is None:
            curr_iter_starting_tree_object = generate_tree_object(starting_tree_path)
        else:
            curr_iter_starting_tree_object = starting_tree_object
        while True:
            curr_iter_run_directory = os.path.join(curr_run_directory , "iter_" + str(spr_iterations_performed_so_far))
            create_or_clean_dir(curr_iter_run_directory)
            curr_iter_starting_tree_path = os.path.join(curr_iter_run_directory , "curr_iter_starting_topology_" + phase_name)
            shutil.copy(curr_iter_starting_tree_topology, curr_iter_starting_tree_path)
            curr_iter_best_topology_path = os.path.join(curr_iter_run_directory , "curr_iter_best_topology_" + phase_name)
            shutil.copy(curr_iter_starting_tree_topology, curr_iter_best_topology_path)
            spr_log_file_object.write("iteration number:" + str(spr_iterations_performed_so_far) + "\n")
            logging.info("iteration number: " + str(spr_iterations_performed_so_far))
            spr_log_file_object.write("running SPR iteration with starting tree:" + curr_iter_starting_tree_path + "\n")
            adjusted_starting_tree_ll=curr_iter_starting_tree_ll
            if phase_name == "use_sampled_MSA":
                adjusted_starting_tree_ll = (curr_iter_starting_tree_ll/INTEGER_CONST) + curr_msa_stats["lasso_intercept"]
            logging.info(
                "iteration raw starting tree ll= {} iteration starting tree adjusted= {} iteration starting tree true ll= {}".format(curr_iter_starting_tree_ll, adjusted_starting_tree_ll,curr_iter_starting_tree_ll_on_full_data))
            curr_iter_best_tree_object, curr_iter_best_ll, curr_iter_best_true_ll, curr_iter_neighbours_tested, curr_iter_ll_comparison_df = SPR_iteration(
                MSA_path,
                curr_msa_stats,
                curr_iter_starting_tree_object,
                curr_iter_starting_tree_ll,
                curr_iter_starting_tree_ll_on_full_data,
                spr_iterations_performed_so_far,
                curr_iter_run_directory,
                curr_iter_best_topology_path,
                spr_log_file_object,
                use_weights)
            if phase_name != "use_sampled_MSA" and abs(curr_iter_best_true_ll - curr_iter_best_ll) > 5:
                logging.error("Problem ! serious different detected between 2 estimations of same full data likelihood")
                print("Problem here!")
            number_of_SPR_neighbours_per_iteration_list.append(curr_iter_neighbours_tested)
            true_vs_sampled_ll_per_iteration_list.append([curr_iter_best_ll, curr_iter_best_true_ll])
            epsilon=EPSILON
            curr_adjusted_ll=curr_iter_best_ll
            if phase_name == "use_sampled_MSA":
                ll_comparison_df = pd.concat([ll_comparison_df, curr_iter_ll_comparison_df])
                epsilon = EPSILON*INTEGER_CONST
                curr_adjusted_ll = (curr_iter_best_ll / INTEGER_CONST) + curr_msa_stats["lasso_intercept"]
            spr_log_file_object.write(
                "iteration summary: log likelihood: {} ; prev iteration ll: {} ; number of SPR neighbours tested: {} \n".format(
                    curr_iter_best_ll,
                    curr_iter_starting_tree_ll, curr_iter_neighbours_tested))
            logging.info(
                "iteration {} raxml raw log likelihood= {} raxml adjusted log likelihood= {} and true log likelihood= {} ".format(
                    spr_iterations_performed_so_far,
                    curr_iter_best_ll, curr_adjusted_ll,
                    curr_iter_best_true_ll))

            if curr_iter_best_ll - curr_iter_starting_tree_ll <= epsilon:
                spr_log_file_object.write(
                    "curr_iteration_ll - prev_iteration_ll <= {}, stopping SPR search\n".format(epsilon))
                logging.info("curr_iteration_ll - prev_iteration_ll <= {}, stopping SPR search\n".format(epsilon))
                break
            ### Updating current iteration results and preparing for next iteration:
            spr_iterations_performed_so_far = spr_iterations_performed_so_far + 1
            curr_iter_starting_tree_topology = curr_iter_best_topology_path
            curr_iter_starting_tree_object = curr_iter_best_tree_object
            curr_iter_starting_tree_ll = curr_iter_best_ll
            curr_iter_starting_tree_ll_on_full_data = curr_iter_best_true_ll
    # delete_file_content(spr_log_file)
    return {
        "current_starting_tree_ll": naive_SPR_search_starting_tree_ll,
        "current_starting_tree_true_ll": naive_SPR_search_starting_tree_true_ll,
        "ll": curr_iter_best_ll,
        "true_ll": curr_iter_best_true_ll,
        "best_topology_newick": curr_iter_best_tree_object.write(format=1),
        "ll_comparison_df": ll_comparison_df,
        "true_vs_sampled_ll_per_iteration_list": true_vs_sampled_ll_per_iteration_list,
        "curr_iter_best_tree_object": curr_iter_best_tree_object,
        "curr_iter_best_topology_path": curr_iter_best_topology_path,
        "spr_moves": spr_iterations_performed_so_far,
        "number_of_SPR_neighbours_per_iteration_list": number_of_SPR_neighbours_per_iteration_list
    }


def get_SPR_file_prefixs(unique_name):
    file_prefix_dict = {}
    file_prefix_dict["run_file_prefix"] = "SPR_"
    file_prefix_dict["ll_on_all_data_run_file_prefix"] = "SPR_"
    return file_prefix_dict


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
        "starting tree is stored in: {}".format(SPR_chosen_starting_tree_path))
    if (full_run):
        logging.info("Starting SPR analysis on full data")
        full_data_param_dict = SPR_search(
            MSA_path=current_file_path,
            run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=SPR_chosen_starting_tree_path,
            starting_tree_object=None,
            phase_name="normal_full_MSA",
            curr_run_directory=curr_run_directory, use_weights = False)
        logging.info("SPR result tree is located in :" + full_data_param_dict["curr_iter_best_topology_path"])
        ll_estimation_sanity_check = abs(full_data_param_dict["ll"] - full_data_param_dict["true_ll"]) > 5
        if ll_estimation_sanity_check:
            logging.error("Likelihood estimation problem detected in full data! ll= {} and true ll= {}  ".format(
                full_data_param_dict["ll"], full_data_param_dict["true_ll"]))
        full_data_SPR_result = {"naive_SPR_ll": full_data_param_dict["ll"],
                                "naive_SPR_spr_moves": full_data_param_dict.get("spr_moves"),
                                "naive_SPR_tree_newick" : full_data_param_dict["best_topology_newick"]
                                # , "full data true vs sampled ll per iteration_list_sanity_check": full_data_param_dict.get(
                                #     "true_vs_sampled_ll_per_iteration_list")
                                }
        logging.info("Full dataset SPR result: " + str(full_data_SPR_result))
        return (full_data_SPR_result)
    else:
        logging.info("Starting SPR analysis on sampled data")
        sub_curr_run_directory = os.path.join(curr_run_directory , "_use_sampled_MSA")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
            #        first_optimized_tree_object, first_optimized_tree_path, ll_comparison_df, true_vs_sampled_ll_per_iteration_list_first_part,
        first_optimized_param_dict = SPR_search(
            MSA_path=curr_msa_stats["sampled_alignment_path"],
            run_unique_name=run_unique_name,
            starting_tree_path=SPR_chosen_starting_tree_path,
            starting_tree_object=None,
            curr_msa_stats=curr_msa_stats,
            curr_run_directory=sub_curr_run_directory,
            phase_name="use_sampled_MSA",
            use_weights = True)
        #use_sampled_starting_tree_ll = first_optimized_param_dict["current_starting_tree_ll"]
        use_sampled_true_starting_tree_ll = first_optimized_param_dict["current_starting_tree_true_ll"]
        tree_newick_first_phase = first_optimized_param_dict["best_topology_newick"]
        ll_comparison_df = first_optimized_param_dict["ll_comparison_df"]
        ll_comparison_df.to_csv(os.path.join(curr_run_directory , "ll_comparison_df.csv"))
        prediction_rho_pearson, prediction_pval_pearson,prediction_rho_spearman, prediction_pval_spearman,mse,mistake_cnt = analyze_ll_comparison_df(ll_comparison_df)
        overall_SPR_steps_list = [first_optimized_param_dict["spr_moves"]]
        first_optimized_tree_path = first_optimized_param_dict["curr_iter_best_topology_path"]
        first_optimized_tree_object = first_optimized_param_dict["curr_iter_best_tree_object"]
        true_vs_sampled_ll_per_iteration_list = [first_optimized_param_dict["true_vs_sampled_ll_per_iteration_list"]]
        overall_SPR_neighbours_per_iteration = [
            first_optimized_param_dict.get("number_of_SPR_neighbours_per_iteration_list")]
        overall_ll_list = [first_optimized_param_dict["true_ll"]]

        ### Continue sampling with full data
        logging.info("Continue SPR analysis using full data")
        # Use previous tree as a starting tree
        full_msa_path = curr_msa_stats["local_alignment_path"]
        sub_curr_run_directory = os.path.join(curr_run_directory ,"_continue_with_full_MSA")
        if not os.path.exists(sub_curr_run_directory):
            os.mkdir(sub_curr_run_directory)
        # SPR_search(MSA_path, run_unique_name, curr_msa_stats, starting_tree_path,starting_tree_object, curr_run_directory, phase_name,
        # use_external_weights=False)
        #  next_best_tree_object, next_optimized_tree, ll_comparison_df_not_usable, true_vs_sampled_ll_per_iteration_list_second_part,

        next_optimized_tree_param_dict = SPR_search(
            MSA_path=full_msa_path, run_unique_name=run_unique_name,
            curr_msa_stats=curr_msa_stats,
            starting_tree_path=first_optimized_tree_path,
            starting_tree_object=first_optimized_tree_object,
            curr_run_directory=sub_curr_run_directory,
            phase_name="continue_with_full_MSA",
            use_weights = False
        )
        ll_estimation_sanity_check = abs(
            next_optimized_tree_param_dict["ll"] - next_optimized_tree_param_dict["true_ll"]) > 5
        if ll_estimation_sanity_check:
            logging.error(
                "Likelihood estimation problem detected in continue with full data! ll= {} and true ll= {}  ".format(
                    next_optimized_tree_param_dict["ll"], next_optimized_tree_param_dict["true_ll"]))
        tree_newick_second_phase =next_optimized_tree_param_dict["best_topology_newick"]
        overall_SPR_steps_list.append(next_optimized_tree_param_dict.get("spr_moves"))
        overall_SPR_neighbours_per_iteration.append(
            next_optimized_tree_param_dict.get("number_of_SPR_neighbours_per_iteration_list"))
        overall_ll_list.append(next_optimized_tree_param_dict.get("ll"))
        true_vs_sampled_ll_per_iteration_list.append(
            next_optimized_tree_param_dict["true_vs_sampled_ll_per_iteration_list"])

        data = {
                "lasso_SPR_first_phase_ll": overall_ll_list[0],
                "lasso_SPR_first_phase_tree_newick" : tree_newick_first_phase,
                "lasso_SPR_first_phase_spr_moves": first_optimized_param_dict["spr_moves"],
                "lasso_SPR_second_phase_ll": overall_ll_list[-1],
                "lasso_SPR_second_phase_tree_newick": tree_newick_second_phase,
                "lasso_SPR_second_phase_spr_moves": next_optimized_tree_param_dict["spr_moves"],
                "R^2_pearson_during_tree_search": prediction_rho_pearson**2,
                "R^2_pearson_during_tree_search_pval": prediction_pval_pearson,
                "spearmanr_during_tree_search": prediction_rho_spearman,
                "spearmanr_during_tree_search_pval": prediction_pval_spearman,
                "mse_during_tree_search" : mse,
                "mistake_cnt" : mistake_cnt,
                "lasso_SPR_starting_tree_ll": use_sampled_true_starting_tree_ll,
                "lasso_SPR_starting_tree_path": SPR_chosen_starting_tree_path

                }
        return data
