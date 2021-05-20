from sklearn.metrics import *
from sklearn import linear_model
from raxml import *
from scipy import stats
import time
import numpy as np





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





def apply_lasso_on_given_data(curr_msa_stats, training_df):
    logging.info("Computing locis and weight using lasso")
    y_training = training_df.sum(axis=1)
    logging.info("Sitelh df dimensions, for lasso computations, are: " + str(training_df.shape))
    start_time = time.time()
    if (curr_msa_stats["alphas"]) != "default":
        if "_" in curr_msa_stats["alphas"]:
            alphas = [float(val) for val in curr_msa_stats["alphas"].split("_")]
        else:
            alphas = [float(curr_msa_stats["alphas"])]
        logging.info("Using given alphas for the Lasso: {alphas}".format(alphas=alphas))
        lasso_model = linear_model.Lasso(normalize=True, max_iter=100000, positive=True, random_state=SEED,
                                           selection='cyclic', alpha=alphas).fit(training_df,
                                                                                  y_training)

    else:
        logging.info("Using default alphas for the Lasso")
        lasso_model = linear_model.LassoCV(cv=2, normalize=True, max_iter=100000, positive=True,
                                           random_state=SEED, selection='cyclic'
                                           ).fit(training_df,
                                                 y_training)

    lasso_training_time = time.time() - start_time
    logging.info("Done training Lasso model. It took {} seconds".format(lasso_training_time))
    return lasso_model, lasso_training_time,y_training





def apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats, curr_run_directory, sitelh_training_df, test_optimized_trees_path):
    step_size = curr_msa_stats["n_partitions"]
    indexes = np.array_split(np.arange(len(sitelh_training_df.columns)), step_size)
    overall_chosen_locis = []
    overall_weights = []
    sampled_alignment_path = os.path.join(curr_run_directory,
                                          curr_msa_stats["file_name"] + "_sampled" + curr_msa_stats[
                                              "file_type"])
    weights_file_path = os.path.join(curr_run_directory, curr_msa_stats.get(
        "file_name") + '_weights.txt')
    lasso_per_batch_csv_path = os.path.join(curr_run_directory, curr_msa_stats.get(
        "file_name") + '_lasso_metrics.csv')
    overall_running_time = 0
    overall_lasso_results = pd.DataFrame()
    for i in range(step_size):
        curr_data = sitelh_training_df.iloc[:,indexes[i]]
        logging.info(f"Applying {i}th batch of Lasso, based on positions {indexes[i]}")
        lasso_model, lasso_training_time, y_training = apply_lasso_on_given_data(curr_msa_stats,curr_data)
        if USE_INTEGER_WEIGHTS:
            weights = [int(lasso_model.coef_[ind] * INTEGER_CONST) for ind in range(len(lasso_model.coef_))]
        else:
            weights = lasso_model.coef_
        chosen_locis = list(np.array([ind for ind in range(len(weights)) if weights[ind] != 0]) + int(indexes[i][0]))
        chosen_loci_weights = [weights[ind] for ind in range(len(weights)) if weights[ind] != 0]
        overall_chosen_locis += chosen_locis
        overall_weights +=chosen_loci_weights
        overall_running_time = overall_running_time+lasso_training_time
        training_r2 = get_training_metrics(lasso_model,curr_data,y_training)[1]["lasso_training_R^2"]
        lasso_results = {
            "training_r2" : training_r2,
                      #"lasso_alpha": lasso_model.alpha_,
                      #"lasso_alphas": lasso_model.alphas_,
                      "n_iter_lasso": lasso_model.n_iter_,
                      "lasso_training_time":lasso_training_time,
                      "lasso_intercept": lasso_model.intercept_,
                      "lasso_training_X": curr_data,
        "lasso_training_y" : y_training,
            "number_loci_chosen" : len(chosen_locis)
        }
        logging.info(f"Results for the {i}th fold: \n {lasso_results}")
        overall_lasso_results = overall_lasso_results.append(lasso_results, ignore_index= True)

    overall_lasso_results.to_csv(lasso_per_batch_csv_path)
    logging.info("Writing overall chosen positions to {}".format(sampled_alignment_path))
    write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], sampled_alignment_path, overall_chosen_locis,
                                    curr_msa_stats["file_type_biopython"])
    logging.info("Writing overall weights to : " + weights_file_path)
    with open(weights_file_path, 'w') as f:
        for weight in overall_weights:
            f.write(str(weight) + " ")
    Lasso_results = ({"number_loci_chosen": len(overall_chosen_locis), "lasso_chosen_locis": overall_chosen_locis,
                      "sample_pct": len(overall_chosen_locis) / curr_msa_stats.get("n_loci"),
                      "lasso_chosen_weights": overall_weights, "weights_file_path": weights_file_path,
                      "sampled_alignment_path": sampled_alignment_path}

    )


    if test_optimized_trees_path is not None:
        test_running_directory = os.path.join(curr_run_directory, "test_ll_evaluations")
        create_dir_if_not_exists(test_running_directory)
        y_test_predicted, y_test_true, test_results = evaluate_lasso_performance_on_test_data(
            test_optimized_trees_path, curr_msa_stats, test_running_directory, sampled_alignment_path,
            weights_file_path, lasso_model.intercept_)
        # if GENERATE_LASSO_DESCRIPTIVE:
        #     generate_lasso_descriptive(y_training_predicted, y_training,
        #                                y_test_predicted, y_test_true, curr_run_directory)
        Lasso_results.update(test_results)
        return Lasso_results




