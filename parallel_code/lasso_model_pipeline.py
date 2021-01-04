
from sklearn import linear_model
from raxml import *
from scipy import stats


def evaluate_lasso_performance_on_test_data(lasso_model, curr_msa_stats):
    logging.info("Evaluating model on test data")
    sitelh_test_df = curr_msa_stats["test_sitelh_df"]
    y_test = sitelh_test_df.sum(axis=1)
    test_predicted_values = lasso_model.predict(sitelh_test_df)
    return test_predicted_values, y_test, sitelh_test_df


def generate_lasso_descriptive(training_predicted_values, y_training, sitelh_training_df,
                               test_predicted_values, y_test, sitelh_test_df,name, curr_run_directory):
    training_sitelh_df_prediction = pd.concat(
        [y_training.rename("y_training"), pd.Series(training_predicted_values).rename("y_training_predicted"),
         sitelh_training_df], axis=1, ignore_index=True, sort=False)
    training_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "training_sitelh_df_prediction_{}.csv".format(name)))
    test_sitelh_df_prediction = pd.concat(
        [y_test.rename("y_test"), pd.Series(test_predicted_values).rename("y_test_predicted"), sitelh_test_df], axis=1,
        ignore_index=True, sort=False)
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction_{}.csv".format(name)))


def apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats,name,curr_run_directory):
    logging.info('Applying Lasso on sitelh data' )
    sitelh_training_df = curr_msa_stats.get("training_sitelh_df")
    lasso_log_file = os.path.join(curr_run_directory,name+"_lasso.log")
    with open(lasso_log_file, 'w') as lasso_log:
            logging.info("Computing locis and weight using lasso")
            y_training = sitelh_training_df.sum(axis=1)
            logging.info("Sitelh df dimensions, for lasso computations, are: " + str(sitelh_training_df.shape))
            lasso_model = linear_model.LassoCV(cv=5, normalize=True, max_iter=100000,positive=True).fit(sitelh_training_df,
                                                                                          y_training)  # add positive=True if using RaxML
            y_training_predicted=lasso_model.predict(sitelh_training_df)
            integer_coefficients = [int(lasso_model.coef_[ind]*INTEGER_CONST) for ind in range(len(lasso_model.coef_))]
            chosen_locis = [ind for ind in range(len(integer_coefficients)) if integer_coefficients[ind] > 0]
            chosen_loci_weights = [integer_coefficients[ind] for ind in range(len(integer_coefficients)) if integer_coefficients[ind] > 0]
            lasso_log.write("Lasso chosen locis are:" + str(chosen_locis) + "\n")
            lasso_log.write("Lasso nonzero coefficient are:" + str(chosen_loci_weights) + "\n")
            training_r_squared=lasso_model.score(sitelh_training_df,y_training)

            training_spearmanr = stats.spearmanr(y_training,y_training_predicted)[0]
            y_test_predicted, y_test, sitelh_test_df = evaluate_lasso_performance_on_test_data(
                lasso_model, curr_msa_stats)
            test_r_squared = lasso_model.score(sitelh_test_df, y_test)
            test_spearmanr = stats.spearmanr(y_test,y_test_predicted)[0]
            if GENERATE_LASSO_DESCRIPTIVE:
                generate_lasso_descriptive(y_training_predicted, y_training, sitelh_training_df,
                                           y_test_predicted, y_test,
                                           sitelh_test_df,name, curr_run_directory)  # lasso_log.write("Lasso R^2 on test data is:" + str(test_r_squared) + "\n")
            weights_file_path = os.path.join(curr_run_directory ,curr_msa_stats.get(
                "file_name") + "_lasso_weights.txt")
            logging.info("About to write weights file to : " + weights_file_path)
            with open(weights_file_path, 'w') as f:
                for weight in chosen_loci_weights:
                    f.write(str(weight) + " ")
            sampled_alignment_path = os.path.join(curr_run_directory,
                                               curr_msa_stats["file_name"] +"_sampled"+ curr_msa_stats["file_type"])
            curr_msa_stats["sampled_alignment_path"] = sampled_alignment_path
            logging.info("Writing only chosen positions to {}".format(sampled_alignment_path))
            write_to_sampled_alingment_path(curr_msa_stats["alignment_data"], sampled_alignment_path, chosen_locis,
                                            curr_msa_stats["file_type_biopython"])
            lambda_function = lambda row: lasso_model.predict(np.reshape(row, (1, -1)))[0]
            samp_indexes_pct=len(chosen_locis) / curr_msa_stats.get("n_loci")
            curr_msa_stats.update({"number_loci_chosen": len(chosen_locis), "lasso_chosen_locis": chosen_locis,"sample_pct": samp_indexes_pct,
                                   "lasso_coeffs": lasso_model.coef_,
                                   "lasso_intercept": lasso_model.intercept_*INTEGER_CONST,
                                   "lasso_chosen_weights": chosen_loci_weights, "weights_file_path": weights_file_path,
                                   "lasso_training_R^2": training_r_squared, "lasso_test_R^2": test_r_squared,
                                   "lasso_training_spearmanr": training_spearmanr,"lasso_test_spearmanr": test_spearmanr,
                                   "lasso_predict_func": lambda_function})
