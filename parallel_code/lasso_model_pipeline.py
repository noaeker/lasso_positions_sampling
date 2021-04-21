from sklearn.metrics import *
from sklearn import linear_model
from raxml import *
from scipy import stats
import pickle
import time
from relaxed_lasso import RelaxedLassoCV


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

    test_r_squared = stats.pearsonr(true_ll_values, lasso_ll_values_adjusted)[0] ** 2
    test_mse = mean_squared_error(true_ll_values, lasso_ll_values_adjusted)
    test_spearmanr = stats.spearmanr(true_ll_values, lasso_ll_values_adjusted)[0]
    results= {"lasso_test_R^2" : test_r_squared, "lasso_test_mse" : test_mse, "lasso_test_spearmanr" : test_spearmanr}
    return lasso_ll_values_adjusted, true_ll_values, results



def generate_lasso_descriptive(training_predicted_values, training_true_values,
                               test_predicted_values, test_true_values, curr_run_directory):
    training_sitelh_df_prediction =   pd.DataFrame({'predicted_training_ll': training_predicted_values ,'true_training_ll': training_true_values})
    training_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "training_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.DataFrame({'predicted_test_ll': test_predicted_values ,'true_test_ll': test_true_values,})
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction.csv"))



def get_training_metrics(lasso_model,sitelh_training_df,y_training):
     y_training_predicted = lasso_model.predict(sitelh_training_df)
     training_r_squared = lasso_model.score(sitelh_training_df, y_training)
     training_mse = mean_squared_error(y_training, y_training_predicted)
     training_spearmanr = stats.spearmanr(y_training, y_training_predicted)[0]
     training_results= {"lasso_training_R^2": training_r_squared,
     "lasso_training_spearmanr": training_spearmanr,
     "lasso_training_mse": training_mse}
     return y_training_predicted, training_results


def extract_required_lasso_results(lasso_model, curr_run_directory, curr_msa_stats, lasso_log):
    if USE_INTEGER_WEIGHTS:
        weights = [int(lasso_model.coef_[ind] * INTEGER_CONST) for ind in range(len(lasso_model.coef_))]
    else:
        weights = lasso_model.coef_
    chosen_locis = [ind for ind in range(len(weights)) if weights[ind] != 0]
    sampled_alignment_path = os.path.join(curr_run_directory,
                                          curr_msa_stats["file_name"] + "_sampled" + curr_msa_stats[
                                              "file_type"])
    logging.info("Writing only chosen positions to {}".format(sampled_alignment_path))
    chosen_loci_weights = [weights[ind] for ind in range(len(weights)) if weights[ind] != 0]
    write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], sampled_alignment_path, chosen_locis,
                                    curr_msa_stats["file_type_biopython"])
    lasso_log.write("Lasso chosen locis are:" + str(chosen_locis) + "\n")
    lasso_log.write("Lasso nonzero coefficient are:" + str(chosen_loci_weights) + "\n")
    lasso_model_file_path = os.path.join(curr_run_directory, "lasso_model.sav")
    pickle.dump(lasso_model, open(lasso_model_file_path, 'wb'))
    weights_file_path = os.path.join(curr_run_directory, curr_msa_stats.get(
        "file_name") + "_lasso_weights.txt")
    logging.info("About to write weights file to : " + weights_file_path)
    with open(weights_file_path, 'w') as f:
        for weight in chosen_loci_weights:
            f.write(str(weight) + " ")
    return chosen_locis, chosen_loci_weights, sampled_alignment_path,lasso_model_file_path,weights_file_path

def apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats, curr_run_directory, sitelh_training_df, test_optimized_trees_path):
    lasso_log_file = os.path.join(curr_run_directory,"lasso.log")
    with open(lasso_log_file, 'w') as lasso_log:
            logging.info("Computing locis and weight using lasso")
            y_training = sitelh_training_df.sum(axis=1)
            logging.info("Sitelh df dimensions, for lasso computations, are: " + str(sitelh_training_df.shape))
            start_time = time.time()
            lasso_model = linear_model.LassoCV(cv=5, normalize=True, max_iter=100000,positive=True, random_state=SEED, selection='cyclic').fit(sitelh_training_df,
                                                                                          y_training)  # add positive=True if using RaxML
            lasso_training_time = time.time()-start_time
            logging.info("Done training Lasso model. It took {} seconds".format(lasso_training_time))
            chosen_locis, chosen_loci_weights, sampled_alignment_path, lasso_model_file_path, weights_file_path = extract_required_lasso_results(lasso_model, curr_run_directory, curr_msa_stats, lasso_log)
            y_training_predicted, training_results = get_training_metrics(lasso_model,sitelh_training_df,y_training)
            test_running_directory = os.path.join(curr_run_directory, "test_ll_evaluations")
            create_dir_if_not_exists(test_running_directory)
            y_test_predicted, y_test_true,test_results = evaluate_lasso_performance_on_test_data(
                test_optimized_trees_path, curr_msa_stats, test_running_directory,sampled_alignment_path,weights_file_path,lasso_model.intercept_)
            if GENERATE_LASSO_DESCRIPTIVE:
                generate_lasso_descriptive(y_training_predicted, y_training,
                                           y_test_predicted, y_test_true, curr_run_directory)

            Lasso_results = ({"number_loci_chosen": len(chosen_locis), "lasso_chosen_locis": chosen_locis,"sample_pct": len(chosen_locis) / curr_msa_stats.get("n_loci"),
                                   "lasso_coeffs": lasso_model.coef_,
                              "lasso_alpha" : lasso_model.alpha_,
                              "n_iter_lasso": lasso_model.n_iter_,
                              "lasso_training_time" :  lasso_training_time,
                                   "lasso_intercept": lasso_model.intercept_,
                                   "lasso_chosen_weights": chosen_loci_weights, "weights_file_path": weights_file_path,
                                   "sampled_alignment_path": sampled_alignment_path,
                                   "lasso_model_file_path":lasso_model_file_path})
            Lasso_results.update(training_results)
            Lasso_results.update(test_results)
            return Lasso_results

