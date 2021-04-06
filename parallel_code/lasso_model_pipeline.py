from sklearn.metrics import *
from sklearn import linear_model
from raxml import *
from scipy import stats
import pickle
import time



def evaluate_lasso_performance_on_test_data(optimized_random_trees_path, curr_msa_stats, curr_run_directory, sampled_alignment_path, weights_file_path, lasso_intercept, opt_brlen = True):
    logging.info("Evaluating model on test optimized random trees")
    local_file_path = curr_msa_stats.get("local_alignment_path")
    prefix_true = "opt_using_full_data"
    true_ll_values =raxml_optimize_trees_for_given_msa(local_file_path, prefix_true, optimized_random_trees_path, curr_msa_stats,
                                                       curr_run_directory, opt_brlen, weights=None, return_trees_file=False)[0]
    prefix_lasso="opt_using_lasso"
    lasso_ll_values = raxml_optimize_trees_for_given_msa(sampled_alignment_path, prefix_lasso, optimized_random_trees_path, curr_msa_stats,
                                                         curr_run_directory, opt_brlen, weights=weights_file_path, return_trees_file=False)[0]
    lasso_ll_values_adjusted = [(ll/INTEGER_CONST)+lasso_intercept for ll in lasso_ll_values]
    return true_ll_values,   lasso_ll_values_adjusted


def generate_lasso_descriptive(training_predicted_values, training_true_values,
                               test_predicted_values, test_true_values, curr_run_directory):
    training_sitelh_df_prediction =   pd.DataFrame({'predicted_training_ll': training_predicted_values ,'true_training_ll': training_true_values})
    training_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "training_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.DataFrame({'predicted_test_ll': test_predicted_values ,'true_test_ll': test_true_values,})
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction.csv"))


def apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats, curr_run_directory, sitelh_training_df, test_optimized_trees_path):
    lasso_log_file = os.path.join(curr_run_directory,"lasso.log")
    with open(lasso_log_file, 'w') as lasso_log:
            logging.info("Computing locis and weight using lasso")
            y_training = sitelh_training_df.sum(axis=1)
            logging.info("Sitelh df dimensions, for lasso computations, are: " + str(sitelh_training_df.shape))
            start_time = time.time()
            lasso_model = linear_model.LassoCV(cv=5, normalize=True, max_iter=100000,positive=True).fit(sitelh_training_df,
                                                                                          y_training)  # add positive=True if using RaxML
            lasso_training_time = time.time()-start_time
            logging.info("Done training Lasso model. It took {} seconds".format(lasso_training_time))
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
                test_optimized_trees_path, curr_msa_stats, test_running_directory,sampled_alignment_path,weights_file_path,lasso_model.intercept_)
            y_test_true_no_brlen_opt, y_test_predicted_no_brlen_opt = evaluate_lasso_performance_on_test_data(
                test_optimized_trees_path, curr_msa_stats, test_running_directory,sampled_alignment_path,weights_file_path,lasso_model.intercept_, opt_brlen= False)
            test_r_squared = stats.pearsonr(y_test_true, y_test_predicted)[0]**2
            test_r_squared_no_brlen_opt = stats.pearsonr(y_test_true_no_brlen_opt, y_test_predicted_no_brlen_opt)[0]**2
            test_mse = mean_squared_error(y_test_true, y_test_predicted)
            test_spearmanr = stats.spearmanr(y_test_true, y_test_predicted)[0]
            if GENERATE_LASSO_DESCRIPTIVE:
                generate_lasso_descriptive(y_training_predicted, y_training,
                                           y_test_predicted, y_test_true,curr_run_directory)

            samp_indexes_pct=len(chosen_locis) / curr_msa_stats.get("n_loci")
            Lasso_results = ({"number_loci_chosen": len(chosen_locis), "lasso_chosen_locis": chosen_locis,"sample_pct": samp_indexes_pct,
                                   "lasso_coeffs": lasso_model.coef_,
                              "lasso_alpha" : lasso_model.alpha_,
                              "n_iter_lasso": lasso_model.n_iter_,
                              "lasso_training_time" :  lasso_training_time,
                                   "lasso_intercept": lasso_model.intercept_,
                                   "lasso_chosen_weights": chosen_loci_weights, "weights_file_path": weights_file_path,
                                   "lasso_training_R^2": training_r_squared, "lasso_test_R^2": test_r_squared,
                                   "lasso_test_R^2_no_brlen_opt": test_r_squared_no_brlen_opt,
                                   "lasso_training_spearmanr": training_spearmanr,"lasso_test_spearmanr": test_spearmanr,
                                   "lasso_training_mse": training_mse,
                                   "lasso_test_mse": test_mse,"sampled_alignment_path": sampled_alignment_path,
                                   "lasso_model_file_path":lasso_model_file_path})
            return Lasso_results

