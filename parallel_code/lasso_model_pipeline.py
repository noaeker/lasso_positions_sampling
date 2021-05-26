from sklearn.metrics import *
from sklearn import linear_model
from raxml import *
from scipy import stats
import time
import numpy as np
from sklearn import preprocessing
import pickle





def evaluate_lasso_performance_on_test_data(curr_msa_stats, curr_run_directory, sampled_alignment_path, weights_file_path, lasso_intercept):
    optimized_random_trees_path = curr_msa_stats["optimized_test_topologies_path"]
    #logging.info("Evaluating model on test optimized random trees")
    true_ll_values = curr_msa_stats["test_ll_values"]
    prefix_lasso="opt_using_lasso"
    lasso_ll_values = raxml_optimize_trees_for_given_msa(sampled_alignment_path, prefix_lasso, optimized_random_trees_path, curr_msa_stats,
                                                         curr_run_directory, opt_brlen=True, weights=weights_file_path, return_trees_file=False)[0]
    lasso_ll_values_adjusted = [(ll/INTEGER_CONST)+lasso_intercept/INTEGER_CONST for ll in lasso_ll_values]
    lasso_ll_values_no_opt = \
    raxml_optimize_trees_for_given_msa(sampled_alignment_path, prefix_lasso, optimized_random_trees_path,
                                       curr_msa_stats,
                                       curr_run_directory, opt_brlen=False, weights=weights_file_path,
                                       return_trees_file=False)[0]
    lasso_ll_values_adjusted_no_opt = [(ll/INTEGER_CONST)+lasso_intercept/INTEGER_CONST for ll in lasso_ll_values_no_opt]

    test_r_squared = stats.pearsonr(true_ll_values, lasso_ll_values_adjusted)[0] ** 2
    test_mse = mean_squared_error(true_ll_values, lasso_ll_values_adjusted)
    test_spearmanr = stats.spearmanr(true_ll_values, lasso_ll_values_adjusted)[0]
    test_r_squared_no_opt = stats.pearsonr(true_ll_values, lasso_ll_values_adjusted_no_opt)[0] ** 2
    test_mse_no_opt = mean_squared_error(true_ll_values, lasso_ll_values_adjusted_no_opt)
    test_spearmanr_no_opt = stats.spearmanr(true_ll_values, lasso_ll_values_adjusted_no_opt)[0]
    results= {"lasso_test_R^2" : test_r_squared, "lasso_test_mse" : test_mse, "lasso_test_spearmanr" : test_spearmanr,
              "lasso_test_R^2_no_opt": test_r_squared_no_opt, "lasso_test_mse_no_opt": test_mse_no_opt, "lasso_test_spearmanr_no_opt": test_spearmanr_no_opt
              }
    return lasso_ll_values_adjusted,lasso_ll_values_adjusted_no_opt, true_ll_values, results



def generate_lasso_descriptive(training_predicted_values,training_predicted_values_no_opt, training_true_values,
                               test_predicted_values, test_true_values, curr_run_directory):
    training_sitelh_df_prediction =   pd.DataFrame({'predicted_training_ll': training_predicted_values ,'true_training_ll': training_true_values})
    training_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "training_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.DataFrame({'predicted_test_ll': test_predicted_values ,'true_test_ll': test_true_values})
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.DataFrame(
        {'predicted_test_ll_no_opt': training_predicted_values_no_opt, 'true_test_ll_no_opt': test_true_values })
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction_no_opt.csv"))



def get_training_metrics(intercept, overall_chosen_locis, overall_weights,sitelh_training_df,y_training):
     y_training_predicted = list((np.array(sitelh_training_df.iloc[:,overall_chosen_locis])).dot((np.array(overall_weights)))+intercept)
     training_r_squared = stats.pearsonr(y_training_predicted, y_training)[0]**2
     training_mse = mean_squared_error(y_training, y_training_predicted)
     training_spearmanr = stats.spearmanr(y_training, y_training_predicted)[0]
     training_results= {"lasso_training_R^2": training_r_squared,
     "lasso_training_spearmanr": training_spearmanr,
     "lasso_training_mse": training_mse}
     return y_training_predicted, training_results




def choose_best_alpha_ind(curr_msa_stats,coeff_path_results):
    for i in range(1,len(coeff_path_results)):
        #curr_test_r2 = coeff_path_results[i]["test_r_2"]
        curr_sample_pct = coeff_path_results[i]["sample_pct"]
        #next_r_2 = [coeff_path_results[j]["test_r_2"] for j in range(i,len(coeff_path_results))]
        if curr_sample_pct>=curr_msa_stats["sample_thresholds"]:
            return i
    return len(coeff_path_results)-1



def generate_weights_file_and_sampled_msa(curr_run_directory, curr_msa_stats,name, chosen_locis,chosen_loci_weights):
    sampled_alignment_path = os.path.join(curr_run_directory,
                                          curr_msa_stats["file_name"] + f"_sampled_{name}" + curr_msa_stats[
                                              "file_type"])
    write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], sampled_alignment_path,
                                    chosen_locis,
                                    curr_msa_stats["file_type_biopython"])
    weights_file_path = os.path.join(curr_run_directory, curr_msa_stats.get(
        "file_name") + f'_weights_{name}.txt')
    with open(weights_file_path, 'w') as f:
        for weight in chosen_loci_weights:
            f.write(str(weight) + " ")
    return sampled_alignment_path, weights_file_path




def evaluate_coeffs_on_test_set(coeffs, ind, alpha, curr_run_directory, curr_msa_stats, y_mean, scaler, calc_r2 = False):
    chosen_locis,chosen_loci_weights = get_chosen_locis_and_weights(coeffs, 0)
    if len(chosen_locis)==0:
        coeff_path_results = {"alpha": alpha, "test_r_2": -1,"sample_pct": 0 }
    else:
        test_running_directory = os.path.join(curr_run_directory, f"test_ll_evaluations_{ind}")
        create_dir_if_not_exists(test_running_directory)
        sampled_alignment_path, weights_file_path = generate_weights_file_and_sampled_msa(test_running_directory, curr_msa_stats,ind, chosen_locis,chosen_loci_weights)
        intercept = y_mean - np.dot(coeffs, scaler.mean_ * np.reciprocal(scaler.scale_))
        coeff_path_results=({"alpha": alpha,
                                   "sample_pct": len(chosen_locis) / curr_msa_stats["n_loci"]})
        if calc_r2:
            test_metrics = evaluate_lasso_performance_on_test_data(
                curr_msa_stats, test_running_directory, sampled_alignment_path,
                weights_file_path, intercept)
            coeff_path_results["test_r_2"] = test_metrics[3]["lasso_test_R^2"]
    return coeff_path_results


def get_best_lasso_model_on_given_data(curr_msa_stats, training_df, curr_run_directory):
    logging.info("Computing locis and weight using lasso")
    y_training = training_df.sum(axis=1)
    logging.info("Sitelh df dimensions, for lasso computations, are: " + str(training_df.shape))
    scaler = preprocessing.StandardScaler().fit(training_df.values)
    training_df_scaled = scaler.transform(training_df)
    y_mean = y_training.mean()
    training_y_scaled = (y_training-y_training.mean())
    lasso_dump_path = os.path.join(curr_run_directory,"lasso_path.dump")
    baseline_lasso_dump_path = lasso_dump_path.replace(curr_msa_stats["run_prefix"],
                                                       curr_msa_stats["lasso_path_baseline_run_prefix"])
    if os.path.exists(baseline_lasso_dump_path):
        logging.info("Using existing lasso path results in {}".format(baseline_lasso_dump_path))
        with open(baseline_lasso_dump_path, 'rb') as handle:
           lasso_model,lasso_training_time = pickle.load(handle)
    else:
        logging.info("Computing Lasso path from beggining")
        start_time = time.time()
        lasso_model = linear_model.lasso_path(X=training_df_scaled,
                                              y=training_y_scaled, eps=1e-7, positive=True, n_alphas=100)
        lasso_training_time = time.time() - start_time
        lasso_dump_data = (lasso_model, lasso_training_time)
        with open(lasso_dump_path, 'wb') as handle:
            pickle.dump(lasso_dump_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    coeffs_path = lasso_model[1]
    alphas = lasso_model[0]
    coeff_path_results = []
    start_time = time.time()
    rescaled_coeffs_path = np.transpose(coeffs_path)*(np.reciprocal(scaler.scale_))
    for ind,coeffs,alpha in zip(range(1,len(alphas)),(rescaled_coeffs_path),alphas[1:]):
        coeffs_metrics = evaluate_coeffs_on_test_set(coeffs, ind, alpha, curr_run_directory, curr_msa_stats, y_mean, scaler)
        coeff_path_results.append( coeffs_metrics)
    test_set_coeff_evaluation_time = time.time()-start_time
    coeff_path_df = pd.DataFrame(coeff_path_results)
    coeff_path_df.to_csv(os.path.join(curr_run_directory,"alphas_vs_r2.csv"))
    chosen_coefficient_ind = choose_best_alpha_ind(curr_msa_stats,coeff_path_results)
    best_coeffs = rescaled_coeffs_path[chosen_coefficient_ind,:]
    alpha = alphas[chosen_coefficient_ind]
    intercept = y_mean- np.dot(coeffs_path [:,chosen_coefficient_ind],scaler.mean_*np.reciprocal(scaler.scale_))
    return best_coeffs ,intercept, lasso_training_time,test_set_coeff_evaluation_time,y_training, alphas, alpha




def get_chosen_locis_and_weights(coeff_array, coef_starting_point):
    if USE_INTEGER_WEIGHTS:
        weights = list(((coeff_array * INTEGER_CONST).astype(int)))
    else:
        weights = list(coeff_array)
    chosen_locis = list(np.array([ind for ind in range(len(weights)) if weights[ind] != 0]) + coef_starting_point)
    chosen_loci_weights = [weights[ind] for ind in range(len(weights)) if weights[ind] != 0]
    return chosen_locis, chosen_loci_weights




def apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats, curr_run_directory, sitelh_training_df, test_optimized_trees_path):
    step_size = curr_msa_stats["n_partitions"]
    indexes = np.array_split(np.arange(len(sitelh_training_df.columns)), step_size)
    overall_chosen_locis = []
    overall_weights = []
    overall_alpha = []
    overall_alphas = []
    overall_lasso_results = pd.DataFrame()
    lasso_per_batch_csv_path = os.path.join(curr_run_directory, curr_msa_stats.get(
        "file_name") + '_lasso_metrics.csv')
    overall_lasso_running_time,overall_test_running_time,overall_intercept  = 0,0,0
    for i in range(step_size):
        curr_data = sitelh_training_df.iloc[:,indexes[i]]
        logging.info(f"Applying {i}th batch of Lasso, based on positions {indexes[i]}")
        chosen_coefficient_array ,intercept, lasso_training_time,test_set_coeff_eval_time,y_training, alphas, alpha = get_best_lasso_model_on_given_data(curr_msa_stats, curr_data, curr_run_directory)
        chosen_locis, chosen_loci_weights = get_chosen_locis_and_weights(chosen_coefficient_array,int(indexes[i][0]))
        if USE_INTEGER_WEIGHTS:
            intercept = intercept*INTEGER_CONST
        overall_intercept = overall_intercept + intercept
        overall_chosen_locis = overall_chosen_locis+chosen_locis
        overall_weights +=chosen_loci_weights
        overall_lasso_running_time = overall_lasso_running_time+lasso_training_time
        overall_test_running_time = overall_test_running_time +test_set_coeff_eval_time
        overall_alpha.append(alpha)
        overall_alphas.append(alphas)
        lasso_results = {
                      "lasso_alpha": overall_alpha,
                      "lasso_alphas": overall_alphas,
                      "lasso_training_time":lasso_training_time,
            "lasso_test_eval_time" : test_set_coeff_eval_time,
                      "lasso_intercept": intercept,
            "number_loci_chosen" : len(chosen_locis)
        }
        logging.info(f"Results for the {i}th fold: \n {lasso_results}")
        overall_lasso_results = overall_lasso_results.append(lasso_results, ignore_index= True)

    overall_lasso_results.to_csv(lasso_per_batch_csv_path)
    sampled_alignment_path,weights_file_path = generate_weights_file_and_sampled_msa(curr_run_directory, curr_msa_stats, "overall", overall_chosen_locis, overall_weights)

    Lasso_results = ({"number_loci_chosen": len(overall_chosen_locis), "lasso_chosen_locis": overall_chosen_locis,
                      "sample_pct": len(overall_chosen_locis) / curr_msa_stats.get("n_loci"),
                      "lasso_chosen_weights": overall_weights, "weights_file_path": weights_file_path,
                      "sampled_alignment_path": sampled_alignment_path, "lasso_intercept" : overall_intercept, "lasso_running_time": overall_lasso_running_time, "lasso_test_eval_time": overall_test_running_time}

    )

    y_training_predicted, training_results = get_training_metrics(overall_intercept, overall_chosen_locis, overall_weights, curr_data, y_training)

    if test_optimized_trees_path is not None:
        test_running_directory = os.path.join(curr_run_directory, "test_ll_evaluations")
        create_dir_if_not_exists(test_running_directory)
        y_test_predicted,y_test_predicted_no_opt, y_test_true, test_results = evaluate_lasso_performance_on_test_data(
           curr_msa_stats, test_running_directory, sampled_alignment_path,
            weights_file_path, overall_intercept)

        if GENERATE_LASSO_DESCRIPTIVE:
             generate_lasso_descriptive(y_training_predicted,y_test_predicted_no_opt, y_training,
                                        y_test_predicted, y_test_true, curr_run_directory)

        Lasso_results.update(test_results)
        print(f'overall test results: {test_results}')
    return Lasso_results




