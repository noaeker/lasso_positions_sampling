from sklearn.metrics import *
from sklearn import linear_model
from raxml import *
from scipy import stats
import pickle



def evaluate_lasso_performance_on_test_data(optimized_random_trees, curr_msa_stats, curr_run_directory,sampled_alignment_path,weights_file_path,lasso_intercept):
    logging.info("Evaluating model on test optimized random trees")
    local_file_path = curr_msa_stats.get("local_alignment_path")
    true_test_ll_vec = []
    lasso_test_ll_vec=[]
    for i in range(len(optimized_random_trees)):
        tree_path=  optimized_random_trees[i]
        prefix_true = "opt_using_full_data_{}".format(i)
        tree_dir = os.path.join(curr_run_directory,"tree_{}".format(i))
        create_dir_if_not_exists(tree_dir)
        true_ll =raxml_optimize_ll_on_given_tree_and_msa(local_file_path, prefix_true, tree_path, curr_msa_stats,
                                                         tree_dir, opt_brlen=True, weights=None, return_tree=False)
        true_test_ll_vec.append(true_ll)
        prefix_lasso="opt_using_lasso_{}".format(i)
        lasso_ll = (raxml_optimize_ll_on_given_tree_and_msa(sampled_alignment_path, prefix_lasso, tree_path, curr_msa_stats,
                                                           tree_dir, opt_brlen=True, weights=weights_file_path, return_tree=False)/INTEGER_CONST)+lasso_intercept
        lasso_test_ll_vec.append(lasso_ll)
    #test_predicted_values_test = [lasso_model.intercept_+(np.dot(chosen_loci_weights,row[chosen_locis]))/INTEGER_CONST for index,row in sitelh_test_df.iterrows()]
    return true_test_ll_vec, lasso_test_ll_vec


def generate_lasso_descriptive(training_predicted_values, training_true_values,
                               test_predicted_values, test_true_values, curr_run_directory):
    training_sitelh_df_prediction =   pd.DataFrame({'predicted_training_ll': training_predicted_values ,'true_training_ll': training_true_values})
    training_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "training_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.DataFrame({'predicted_test_ll': test_predicted_values ,'true_test_ll': test_true_values,})
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction.csv"))


def apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats, curr_run_directory, sitelh_training_df, test_optimized_trees):
    logging.info('Applying Lasso on sitelh data' )
    lasso_log_file = os.path.join(curr_run_directory,"lasso.log")
    with open(lasso_log_file, 'w') as lasso_log:
            logging.info("Computing locis and weight using lasso")
            y_training = sitelh_training_df.sum(axis=1)
            logging.info("Sitelh df dimensions, for lasso computations, are: " + str(sitelh_training_df.shape))
            lasso_model = linear_model.LassoCV(cv=5, normalize=True, max_iter=100000,positive=True).fit(sitelh_training_df,
                                                                                          y_training)  # add positive=True if using RaxML
            lasso_model_file_path = os.path.join(curr_run_directory,"lasso_model.sav")
            pickle.dump(lasso_model, open(lasso_model_file_path, 'wb'))
            y_training_predicted=lasso_model.predict(sitelh_training_df)
            if USE_INTEGER_WEIGHTS:
                weights = [int(lasso_model.coef_[ind]*INTEGER_CONST) for ind in range(len(lasso_model.coef_))]
            else:
                weights = lasso_model.coef_
            chosen_locis = [ind for ind in range(len(weights)) if weights[ind] !=0]
            chosen_loci_weights = [weights[ind] for ind in range(len(weights)) if weights[ind] !=0]
            lasso_log.write("Lasso chosen locis are:" + str(chosen_locis) + "\n")
            lasso_log.write("Lasso nonzero coefficient are:" + str(chosen_loci_weights) + "\n")
            training_r_squared=lasso_model.score(sitelh_training_df,y_training)
            training_mse = mean_squared_error(y_training, y_training_predicted)
            training_spearmanr = stats.spearmanr(y_training,y_training_predicted)[0]
            weights_file_path = os.path.join(curr_run_directory ,curr_msa_stats.get(
            "file_name") + "_lasso_weights.txt")
            logging.info("About to write weights file to : " + weights_file_path)
            with open(weights_file_path, 'w') as f:
                for weight in chosen_loci_weights:
                    f.write(str(weight) + " ")
            sampled_alignment_path = os.path.join(curr_run_directory,
                                               curr_msa_stats["file_name"] +"_sampled"+ curr_msa_stats["file_type"])
            logging.info("Writing only chosen positions to {}".format(sampled_alignment_path))
            write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], sampled_alignment_path, chosen_locis,
                                            curr_msa_stats["file_type_biopython"])
            test_running_directory = os.path.join(curr_run_directory, "test_ll_evaluations")
            create_dir_if_not_exists(test_running_directory)
            y_test_true, y_test_predicted = evaluate_lasso_performance_on_test_data(
                test_optimized_trees, curr_msa_stats, test_running_directory,sampled_alignment_path,weights_file_path,lasso_model.intercept_)
            test_r_squared = stats.pearsonr(y_test_true, y_test_predicted)[0]**2
            test_mse = mean_squared_error(y_test_true, y_test_predicted)
            test_spearmanr = stats.spearmanr(y_test_true, y_test_predicted)[0]
            if GENERATE_LASSO_DESCRIPTIVE:
                generate_lasso_descriptive(y_training_predicted, y_training,
                                           y_test_predicted, y_test_true,curr_run_directory)

            samp_indexes_pct=len(chosen_locis) / curr_msa_stats.get("n_loci")
            Lasso_results = ({"number_loci_chosen": len(chosen_locis), "lasso_chosen_locis": chosen_locis,"sample_pct": samp_indexes_pct,
                                   "lasso_coeffs": lasso_model.coef_,
                                   "lasso_intercept": lasso_model.intercept_,
                                   "lasso_chosen_weights": chosen_loci_weights, "weights_file_path": weights_file_path,
                                   "lasso_training_R^2": training_r_squared, "lasso_test_R^2": test_r_squared,
                                   "lasso_training_spearmanr": training_spearmanr,"lasso_test_spearmanr": test_spearmanr,
                                   "lasso_training_mse": training_mse,
                                   "lasso_test_mse": test_mse,"sampled_alignment_path": sampled_alignment_path,
                                   "lasso_model_file_path":lasso_model_file_path})
            return Lasso_results

