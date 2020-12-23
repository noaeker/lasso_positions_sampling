
from sklearn import linear_model
from code.raxml import *


def evaluate_lasso_performance_on_test_data(lasso_model, curr_msa_stats, n_iter):
    logging.info("Evaluating model on test data")
    sitelh_test_ll_list = []
    for i in range(n_iter):
        curr_site_lh = raxml_compute_per_site_ll_on_a_random_tree(curr_msa_stats["local_alignment_path"], i,
                                                                  curr_msa_stats)
        sitelh_test_ll_list.append(curr_site_lh)
    sitelh_test_df = pd.DataFrame(sitelh_test_ll_list, columns=list(range(len(sitelh_test_ll_list[0]))),
                                  index=list(range(len(sitelh_test_ll_list))))
    y_test = sitelh_test_df.sum(axis=1)
    test_r_squared = lasso_model.score(sitelh_test_df, y_test)
    test_predicted_values = lasso_model.predict(sitelh_test_df)
    return test_r_squared, test_predicted_values, y_test, sitelh_test_df


def generate_lasso_descriptive(curr_msa_stats, training_predicted_values, y_training, sitelh_training_df,
                               test_predicted_values, y_test, sitelh_test_df):
    training_sitelh_df_prediction = pd.concat(
        [y_training.rename("y_training"), pd.Series(training_predicted_values).rename("y_training_predicted"),
         sitelh_training_df], axis=1, ignore_index=True, sort=False)
    training_sitelh_df_prediction.to_csv(
        os.path.join(curr_msa_stats["curr_msa_version_folder"], "training_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.concat(
        [y_test.rename("y_test"), pd.Series(test_predicted_values).rename("y_test_predicted"), sitelh_test_df], axis=1,
        ignore_index=True, sort=False)
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_msa_stats["curr_msa_version_folder"], "test_sitelh_df_prediction.csv"))


def apply_lasso_on_data_and_update_stats(curr_msa_stats):
    logging.info('Applying Lasso on sitelh data' )
    sitelh_training_df = curr_msa_stats.get("sitelh_df")
    lasso_log_file = os.path.join(curr_msa_stats["curr_msa_version_folder"] , "lasso.log")
    with open(lasso_log_file, 'w') as lasso_log:
            logging.info("Computing locis and weight using lasso")
            y_training = sitelh_training_df.sum(axis=1)
            logging.info("Sitelh df dimensions, for lasso computations, are: " + str(sitelh_training_df.shape))
            lasso_model = linear_model.LassoCV(cv=5, normalize=True, max_iter=100000).fit(sitelh_training_df,
                                                                                          y_training)  # add positive=True if using RaxML
            chosen_locis = [ind for ind in range(len(lasso_model.coef_)) if abs(lasso_model.coef_[ind]) > 0]
            chosen_loci_weights = [(lasso_model.coef_[ind]) for ind in chosen_locis]
            lasso_log.write("Lasso chosen locis are:" + str(chosen_locis) + "\n")
            lasso_log.write("Lasso nonzero coefficient are:" + str(chosen_loci_weights) + "\n")
            training_r_squared, training_predicted_values = lasso_model.score(sitelh_training_df,
                                                                              y_training), lasso_model.predict(
                sitelh_training_df)
            test_r_squared, test_predicted_values, y_test, sitelh_test_df = evaluate_lasso_performance_on_test_data(
                lasso_model, curr_msa_stats, RANDOM_TREES_TEST_SIZE)
            if GENERATE_LASSO_DESCRIPTIVE:
                generate_lasso_descriptive(curr_msa_stats, training_predicted_values, y_training, sitelh_training_df,
                                           test_predicted_values, y_test,
                                           sitelh_test_df)  # lasso_log.write("Lasso R^2 on test data is:" + str(test_r_squared) + "\n")
            weights_file_path = os.path.join(curr_msa_stats.get("curr_msa_version_folder") ,curr_msa_stats.get(
                "file_name") + "_lasso_weights.txt")
            logging.info("About to write weights file to : " + weights_file_path)
            with open(weights_file_path, 'w') as f:
                for weight in chosen_loci_weights:
                    f.write(str(weight) + " ")
            lambda_function = lambda row: lasso_model.predict(np.reshape(row, (1, -1)))[0]
            samp_indexes_pct=len(chosen_locis) / curr_msa_stats.get("n_loci")
            curr_msa_stats.update({"number_loci_chosen": len(chosen_locis), "lasso_chosen_locis": chosen_locis,"sample_pct": samp_indexes_pct,
                                   "lasso_coeffs": lasso_model.coef_,
                                   "lasso_intercept": lasso_model.intercept_,
                                   "lasso_chosen_weights": chosen_loci_weights, "weights_file_path": weights_file_path,
                                   "lasso_training_R^2": training_r_squared, "lasso_test_R^2": test_r_squared,
                                   "lasso_predict_func": lambda_function})
