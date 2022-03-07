from sklearn.metrics import *
from sklearn import linear_model
from raxml import *
from scipy import stats
import time
import numpy as np
from sklearn import preprocessing
import itertools
import pickle



def evaluate_lasso_performance_on_test_data(curr_msa_stats, curr_run_directory, sampled_alignment_path,
                                            weights_file_path, lasso_intercept,n_chosen_locis, test_data = None, partitioned_model = None):
    if not test_data:
        test_data = curr_msa_stats
    optimized_random_trees_path = test_data.get("optimized_test_topologies_path",test_data.get("test_optimized_trees_path"))
    if curr_msa_stats.get("pars_optimized_model"):
            edit_num_locis_in_model_file_no_partition(curr_msa_stats["pars_optimized_model"], n_chosen_locis)

    # logging.info("Evaluating model on test optimized random trees")
    true_ll_values = test_data["test_ll_values"]
    prefix_lasso = "opt_using_lasso"
    lasso_ll_values = \
        raxml_optimize_trees_for_given_msa(sampled_alignment_path, prefix_lasso, optimized_random_trees_path,
                                           curr_msa_stats,
                                           curr_run_directory, opt_brlen=True, weights=weights_file_path,
                                           return_trees_file=False)[0]
    lasso_ll_values_adjusted = [(ll / INTEGER_CONST) + lasso_intercept / INTEGER_CONST for ll in lasso_ll_values]
    prefix_lasso = "eval_using_lasso"
    lasso_ll_values_no_opt = \
        raxml_optimize_trees_for_given_msa(sampled_alignment_path, prefix_lasso, optimized_random_trees_path,
                                           curr_msa_stats,
                                           curr_run_directory, opt_brlen=False, weights=weights_file_path,
                                           return_trees_file=False)[0]
    lasso_ll_values_adjusted_no_opt = [(ll / INTEGER_CONST) + lasso_intercept / INTEGER_CONST for ll in
                                       lasso_ll_values_no_opt]

    test_r_squared = stats.pearsonr(true_ll_values, lasso_ll_values_adjusted)[0] ** 2
    test_mse = mean_squared_error(true_ll_values, lasso_ll_values_adjusted)
    test_spearmanr = stats.spearmanr(true_ll_values, lasso_ll_values_adjusted)[0]
    test_r_squared_no_opt = stats.pearsonr(true_ll_values, lasso_ll_values_adjusted_no_opt)[0] ** 2
    test_mse_no_opt = mean_squared_error(true_ll_values, lasso_ll_values_adjusted_no_opt)
    test_spearmanr_no_opt = stats.spearmanr(true_ll_values, lasso_ll_values_adjusted_no_opt)[0]
    results = {"lasso_test_R^2": test_r_squared, "lasso_test_mse": test_mse, "lasso_test_spearmanr": test_spearmanr,
               "lasso_test_R^2_no_opt": test_r_squared_no_opt, "lasso_test_mse_no_opt": test_mse_no_opt,
               "lasso_test_spearmanr_no_opt": test_spearmanr_no_opt
               }
    return lasso_ll_values_adjusted, lasso_ll_values_adjusted_no_opt, true_ll_values, results


def generate_lasso_descriptive(training_predicted_values, training_true_values,
                               test_predicted_values, test_predicted_values_no_opt, test_true_values,
                               curr_run_directory):
    training_sitelh_df_prediction = pd.DataFrame(
        {'predicted_training_ll': training_predicted_values, 'true_training_ll': training_true_values})
    training_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "training_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.DataFrame(
        {'predicted_test_ll': test_predicted_values, 'true_test_ll': test_true_values})
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction.csv"))
    test_sitelh_df_prediction = pd.DataFrame(
        {'predicted_test_ll_no_opt': test_predicted_values_no_opt, 'true_test_ll_no_opt': test_true_values})
    test_sitelh_df_prediction.to_csv(
        os.path.join(curr_run_directory, "test_sitelh_df_prediction_no_opt.csv"))


def get_training_metrics(intercept, overall_chosen_locis, overall_weights, sitelh_training_df, y_training):
    y_training_predicted = list(
        (np.array(sitelh_training_df.iloc[:, overall_chosen_locis])).dot((np.array(overall_weights))) + intercept)
    training_r_squared = stats.pearsonr(y_training_predicted, y_training)[0] ** 2
    training_mse = mean_squared_error(y_training, y_training_predicted)
    training_spearmanr = stats.spearmanr(y_training, y_training_predicted)[0]
    training_results = {"lasso_training_R^2": training_r_squared,
                        "lasso_training_spearmanr": training_spearmanr,
                        "lasso_training_mse": training_mse}
    return y_training_predicted, training_results


def choose_coeffs_ind_for_given_threshold(coeff_path_results, threshold):
    for i in range(1, len(coeff_path_results)):
        # curr_test_r2 = coeff_path_results[i]["test_r_2"]
        curr_sample_pct = coeff_path_results[i]["sample_pct"]
        # next_r_2 = [coeff_path_results[j]["test_r_2"] for j in range(i,len(coeff_path_results))]
        if curr_sample_pct >= threshold:
            return i
    return None


def generate_weights_file_and_sampled_msa(curr_run_directory, curr_msa_stats, name, chosen_locis, chosen_loci_weights):
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


def evaluate_coeffs_on_test_set(coeffs, ind, lambd, curr_run_directory, curr_msa_stats, y_mean, scaler, calc_r2=False):
    chosen_locis, chosen_loci_weights = get_chosen_locis_and_weights(coeffs, 0)
    if len(chosen_locis) == 0:
        coeff_path_results = {"lambd": lambd, "test_r_2": -1, "sample_pct": 0}
    else:
        test_running_directory = os.path.join(curr_run_directory, f"test_ll_evaluations_{ind}")
        create_dir_if_not_exists(test_running_directory)
        sampled_alignment_path, weights_file_path = generate_weights_file_and_sampled_msa(test_running_directory,
                                                                                          curr_msa_stats, ind,
                                                                                          chosen_locis,
                                                                                          chosen_loci_weights)
        intercept = y_mean - np.dot(coeffs, scaler.mean_ * np.reciprocal(scaler.scale_))
        coeff_path_results = ({"lambd": lambd,
                               "sample_pct": len(chosen_locis) / curr_msa_stats["n_loci"]})
        if calc_r2:
            test_metrics = evaluate_lasso_performance_on_test_data(
                curr_msa_stats, test_running_directory, sampled_alignment_path,
                weights_file_path, intercept, len(chosen_locis))
            coeff_path_results["test_r_2"] = test_metrics[3]["lasso_test_R^2"]
    return coeff_path_results


def get_sklearn_lasso_path_on_given_data(curr_msa_stats, training_df, curr_run_directory):
    y_training = training_df.sum(axis=1)
    logging.debug("   ***Sitelh df dimensions, for lasso computations, are: " + str(training_df.shape))
    scaler = preprocessing.StandardScaler().fit(training_df.values)
    training_df_scaled = scaler.transform(training_df)
    y_mean = y_training.mean()
    training_y_scaled = (y_training - y_training.mean())
    start_time = time.time()
    selection = 'random' if curr_msa_stats["random_lasso"] else 'cyclic'
    lasso_model = linear_model.lasso_path(X=training_df_scaled,
                                          y=training_y_scaled, eps=1e-7, positive=True, n_lambds=100,
                                          selection=selection, random_state=SEED)
    lasso_training_time = time.time() - start_time
    coeffs_path = lasso_model[1]
    lambds = lasso_model[0]
    coeff_path_results = []
    start_time = time.time()
    rescaled_coeffs_path = np.transpose(coeffs_path) * (np.reciprocal(scaler.scale_))
    lasso_path_metrics_folder = os.path.join(curr_run_directory, f"lasso_path_metrics_folder")
    create_dir_if_not_exists(lasso_path_metrics_folder)
    for ind, coeffs, lambd in zip(range(1, len(lambds)), (rescaled_coeffs_path), lambds[1:]):
        coeffs_metrics = evaluate_coeffs_on_test_set(coeffs, ind, lambd, lasso_path_metrics_folder, curr_msa_stats,
                                                     y_mean,
                                                     scaler)
        coeff_path_results.append(coeffs_metrics)
    test_set_coeff_evaluation_time = time.time() - start_time
    coeff_path_df = pd.DataFrame(coeff_path_results)
    coeff_path_df.to_csv(os.path.join(curr_run_directory, "lambds_vs_r2.csv"))
    return {"coeff_path_results": coeff_path_results, "lambds": lambds, "rescaled_coeffs_path": rescaled_coeffs_path,
            "coeffs_path": coeffs_path, "y_mean": y_mean, "scaler": scaler, "lasso_training_time": lasso_training_time,
            "test_set_coeff_evaluation_time": test_set_coeff_evaluation_time, "y_training": y_training}
    # return best_coeffs ,intercept, lasso_training_time,test_set_coeff_evaluation_time,y_training, lambds, lambd


def get_sklearn_coeffs_for_given_threshold(threshold, lasso_path_results):
    chosen_coefficient_ind = choose_coeffs_ind_for_given_threshold(lasso_path_results["coeff_path_results"], threshold)
    if chosen_coefficient_ind is None:
        return None, None, None
    best_coeffs = lasso_path_results["rescaled_coeffs_path"][chosen_coefficient_ind, :]
    lambd = lasso_path_results["lambds"][chosen_coefficient_ind]
    intercept = lasso_path_results["y_mean"] - np.dot(lasso_path_results["coeffs_path"][:, chosen_coefficient_ind],
                                                      lasso_path_results["scaler"].mean_ * np.reciprocal(
                                                          lasso_path_results["scaler"].scale_))
    return best_coeffs, intercept, lambd


def get_glmnet_coeffs_for_given_threshold(threshold, MSA_n_loci, glmnet_lasso_path):
    relevant_data = glmnet_lasso_path[glmnet_lasso_path["desired_pcts"] == threshold]
    if relevant_data.empty:
        return None, None, None
    t_intercept = float(relevant_data[relevant_data["loci"] == "(Intercept)"].iloc[0]["estimate"])
    t_lambda = max(relevant_data["lambda"])
    t_coefficients_data = relevant_data[relevant_data["loci"] != "(Intercept)"][["loci", "estimate"]]
    t_coefficients = np.zeros(MSA_n_loci)
    np.put(t_coefficients, t_coefficients_data["loci"], t_coefficients_data["estimate"])
    return t_coefficients, t_intercept, t_lambda


def get_chosen_locis_and_weights(coeff_array, coef_starting_point):
    if USE_INTEGER_WEIGHTS:
        weights = list(((coeff_array * INTEGER_CONST).astype(int)))
    else:
        weights = list(coeff_array)
    chosen_locis = list(np.array([ind for ind in range(len(weights)) if weights[ind] != 0]) + coef_starting_point)
    chosen_loci_weights = [weights[ind] for ind in range(len(weights)) if weights[ind] != 0]
    return chosen_locis, chosen_loci_weights


def unify_msa_and_weights(results_df_per_threshold_and_partition, curr_run_directory, curr_msa_stats,
                          sitelh_training_df, test_optimized_trees_path):
    outputs_per_threshold = {}
    y_training = sitelh_training_df.sum(axis=1)
    logging.debug("Unifying MSAs and weights for each partition")
    for threshold in list(results_df_per_threshold_and_partition["threshold"]):
        threshold_folder = os.path.join(curr_run_directory, f"threshold_{threshold}_outputs")
        create_dir_if_not_exists(threshold_folder)
        t_data = results_df_per_threshold_and_partition[
            results_df_per_threshold_and_partition["threshold"] == threshold]
        t_data.sort_values("partition", inplace=True)
        t_weights = list(itertools.chain.from_iterable(list(t_data["chosen_weights"])))
        t_chosen_locis = list(itertools.chain.from_iterable(list(t_data["chosen_locis"])))
        t_intercept = sum(t_data["lasso_intercept"])
        t_running_time = sum(t_data["lasso_running_time"])
        t_lambds = list(t_data["lambd"])
        t_sampled_alignment_path, t_weights_file_path = generate_weights_file_and_sampled_msa(threshold_folder,
                                                                                              curr_msa_stats,
                                                                                              f't_{threshold}',
                                                                                              t_chosen_locis,
                                                                                              t_weights)
        t_lasso_results = ({"number_loci_chosen": len(t_chosen_locis), "lasso_chosen_locis": t_chosen_locis,
                            "sample_pct": len(t_chosen_locis) / curr_msa_stats.get("n_loci"),
                            "lasso_chosen_weights": t_weights, "weights_file_path": t_weights_file_path,
                            "lambds": t_lambds,
                            "sampled_alignment_path": t_sampled_alignment_path, "lasso_intercept": t_intercept,
                            "lasso_running_time": t_running_time, "lasso_threshold": threshold}

        )

        with open(t_lasso_results["sampled_alignment_path"]) as sampled_path:
            sampled_data = list(SeqIO.parse(sampled_path, curr_msa_stats["file_type_biopython"]))
        sampled_alignment_df = alignment_list_to_df(sampled_data)
        constant_sites_pct, avg_entropy, gap_positions_pct = get_positions_stats(sampled_alignment_df)
        try:
            rate4_site_values = [curr_msa_stats["rate4site_scores"][i] for i in t_lasso_results["lasso_chosen_locis"]]
        except Exception as e:
            logging.error('Failed to get Lasso rate4site scores: ' + str(e))
            rate4_site_values = [-1]
        mean_rate4_site = np.mean(rate4_site_values)
        results_dict = {"lasso_constant_sites_pct": constant_sites_pct, "lasso_avg_entropy": avg_entropy,
                        "lasso_gap_positions_pct": gap_positions_pct, "lasso_rates_4_site": rate4_site_values, "lasso"
                                                                                                               "_mean_rate4site": mean_rate4_site}
        t_lasso_results.update(results_dict)
        lasso_rate4_site_values_path = os.path.join(curr_run_directory, "lasso_rate_4_site")
        try:
            with open(lasso_rate4_site_values_path, 'w') as LASSO_RATE4SITE:
                LASSO_RATE4SITE.write(" ".join([str(r) for r in rate4_site_values]))
        except Exception as e:
            logging.error('Failed to write Lasso rate4site scores: ' + str(e))
        y_training_predicted, training_results = get_training_metrics(t_intercept, t_chosen_locis, t_weights,
                                                                      sitelh_training_df, y_training)

        if test_optimized_trees_path is not None:
            test_running_directory = os.path.join(curr_run_directory, "test_ll_evaluations")
            create_dir_if_not_exists(test_running_directory)
            y_test_predicted, y_test_predicted_no_opt, y_test_true, test_results = evaluate_lasso_performance_on_test_data(
                curr_msa_stats, test_running_directory, t_sampled_alignment_path,
                t_weights_file_path, t_intercept, len(t_chosen_locis))

            if GENERATE_LASSO_DESCRIPTIVE:
                generate_lasso_descriptive(y_training_predicted, y_training,
                                           y_test_predicted, y_test_predicted_no_opt, y_test_true, threshold_folder)
            t_lasso_results.update(test_results)
            t_lasso_results_print = {key: t_lasso_results[key] for key in t_lasso_results if
                                     key not in IGNORE_COLS_IN_CSV}

            logging.info(f"   ***Unified results for threshold : {threshold} are: \n {t_lasso_results_print}")
        outputs_per_threshold[threshold] = t_lasso_results
    return outputs_per_threshold

    return outputs_per_threshold


# generate_lasso_descriptive(training_predicted_values, training_predicted_values_no_opt, training_true_values,
#                               test_predicted_values, test_true_values, curr_run_directory)


def get_glmnet_lasso_path_on_given_data(curr_msa_stats, curr_data, partition_folder, i):
    relaxed_lasso = 1 if curr_msa_stats["relaxed_lasso"] else 0
    curr_data_path = os.path.join(partition_folder, f"partition_{i}_sitelh.csv")
    curr_data.to_csv(curr_data_path, index=False)
    lasso_thresholds = curr_msa_stats['lasso_thresholds']
    command = f"module load R/3.6.1;Rscript --vanila {R_CODE_PATH} {curr_data_path} {partition_folder} {relaxed_lasso} {lasso_thresholds}"
    logging.info(f"About to run lasso command in glmnet: {command}")
    lasso_start_time = time.time()
    # os.system('module load R/3.6.1')
    os.system(command)
    logging.info("R glmnet command is done!")

    lasso_output_file_path = os.path.join(partition_folder, "r_lasso_relaxed.csv") if curr_msa_stats[
        "relaxed_lasso"] else os.path.join(partition_folder, "r_lasso.csv")
    glmnet_lasso_path = pd.read_csv(lasso_output_file_path)
    glmnet_running_time = time.time() - lasso_start_time
    logging.info(f"Lasso results should be found in : {lasso_output_file_path} ")
    return glmnet_lasso_path, glmnet_running_time


def apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats, curr_run_directory, sitelh_training_df,
                                                     test_optimized_trees_path):
    n_partitions = curr_msa_stats["n_partitions"]
    partition_indexes = np.array_split(np.arange(len(sitelh_training_df.columns)), n_partitions)
    lasso_results_per_partition_and_threshold = pd.DataFrame()
    for i in range(n_partitions):
        partition_folder = os.path.join(curr_run_directory, f"partition_{i}_results")
        create_dir_if_not_exists(partition_folder)
        curr_data = sitelh_training_df.iloc[:, partition_indexes[i]]
        logging.debug(f"Applying {i}th batch of Lasso, based on positions {partition_indexes[i]}")
        if curr_msa_stats["use_glmnet_lasso"]:
            glmnet_lasso_path, glmnet_lasso_running_time = get_glmnet_lasso_path_on_given_data(curr_msa_stats,
                                                                                               curr_data,
                                                                                               partition_folder, i)
        else:
            sklearn_lasso_path = get_sklearn_lasso_path_on_given_data(curr_msa_stats, curr_data, partition_folder)
        lasso_thresholds = [float(t) for t in curr_msa_stats['lasso_thresholds'].split("_")]
        for threshold in lasso_thresholds:
            if curr_msa_stats["use_glmnet_lasso"]:
                t_coefficients, t_intercept, t_lambda = get_glmnet_coeffs_for_given_threshold(threshold,
                                                                                              curr_msa_stats["n_loci"],
                                                                                              glmnet_lasso_path)
            else:
                t_coefficients, t_intercept, t_lambda = get_sklearn_coeffs_for_given_threshold(threshold,
                                                                                               sklearn_lasso_path)
            if t_coefficients is None:
                continue
            t_chosen_locis, t_chosen_loci_weights = get_chosen_locis_and_weights(t_coefficients,
                                                                                 int(partition_indexes[i][0]))

            if USE_INTEGER_WEIGHTS:
                t_intercept = t_intercept * INTEGER_CONST
            t_lasso_metrics = {"threshold": threshold,
                               "partition": i,
                               "chosen_locis": t_chosen_locis,
                               "chosen_weights": t_chosen_loci_weights,
                               "lambd": t_lambda,
                               "lasso_intercept": t_intercept,
                               "number_loci_chosen": len(t_chosen_locis),
                               "lasso_running_time": glmnet_lasso_running_time if curr_msa_stats[
                                   "use_glmnet_lasso"] else sklearn_lasso_path["lasso_training_time"],

                               }
            t_lasso_metrics_print = {key: t_lasso_metrics[key] for key in
                                     ["threshold", "partition", "number_loci_chosen", "lambd", "lasso_intercept"]}
            logging.debug(f"Results for the {i}th fold on threshold {threshold}: \n {t_lasso_metrics_print}")
            lasso_results_per_partition_and_threshold = lasso_results_per_partition_and_threshold.append(
                t_lasso_metrics, ignore_index=True)
    lasso_results_per_partition_and_threshold.to_csv(
        os.path.join(curr_run_directory, "lasso_results_per_partition_and_threshold.csv"), index=False)

    outputs_per_threshold = unify_msa_and_weights(lasso_results_per_partition_and_threshold, curr_run_directory,
                                                  curr_msa_stats, sitelh_training_df, test_optimized_trees_path)

    return outputs_per_threshold


def compare_lasso_to_naive_approaches_on_test_set(curr_msa_stats, curr_run_directory, threshold,n_random_iterations = 5):
    comparisons_folder = os.path.join(curr_run_directory,"comparisons")
    os.mkdir(comparisons_folder)
    random_results = pd.DataFrame()
    test_dump_path = os.path.join(curr_run_directory,f"Lasso_folder/test_{curr_msa_stats['random_trees_test_size']}_random_trees_eval/test_set.dump")
    if not os.path.exists(test_dump_path):
        test_dump_path = test_dump_path.replace(curr_msa_stats["run_prefix"],
                                                   curr_msa_stats["training_set_baseline_run_prefix"])

    test_data = pickle.load(open(test_dump_path, "rb"))
    for i in range(n_random_iterations):
        output_dict = {}
        sampled_alignment_path = os.path.join(comparisons_folder,f"random_{i}")
        n_loci = curr_msa_stats["actual_n_loci"]
        random.seed(SEED+i)
        samp_indexes = np.random.choice(range(n_loci), size = int(threshold*n_loci))
        write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], sampled_alignment_path, samp_indexes, file_type = 'fasta')
        curr_random_results = evaluate_lasso_performance_on_test_data(curr_msa_stats, comparisons_folder, sampled_alignment_path,
                                                weights_file_path = None, lasso_intercept = 0,test_data = test_data, n_chosen_locis = int(threshold*n_loci))[3]
        random_results = random_results.append(curr_random_results, ignore_index= True)
    mean_random_results = random_results.mean().to_dict()
    update_dict_with_a_suffix(output_dict, mean_random_results, suffix="_random")
    rates = curr_msa_stats["rate4site_scores"]
    locis = list(range(n_loci))
    locis.sort(key = lambda x: rates[x], reverse= True)
    hightest_rates = locis[:int(threshold*n_loci)]
    sampled_alignment_path = os.path.join(comparisons_folder, f"highest_Rates_{i}")
    write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], sampled_alignment_path, hightest_rates,
                                    file_type='fasta')
    high_rate_results = evaluate_lasso_performance_on_test_data(curr_msa_stats, comparisons_folder, sampled_alignment_path,
                                                             weights_file_path=None, lasso_intercept=0, n_chosen_locis = int(threshold*n_loci), test_data = test_data)[3]
    update_dict_with_a_suffix(output_dict,high_rate_results, suffix="_high_rate")
    shutil.rmtree(comparisons_folder)
    return output_dict


