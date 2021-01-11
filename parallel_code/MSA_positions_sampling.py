from lasso_model_pipeline import *
from generate_SPR import *
from training_and_test_set_generation import *
from raxml import *


def basic_pipeline_for_curr_starting_tree(curr_msa_stats, i, starting_tree_type,
                                          tree_folder, starting_tree_path):
    logging.info("About to start SPR pipeline on random tree {i}".format(i=i))
    logging.info("Current starting tree folder:{folder}\nCurrent starting tree type:{type}\nCurrent starting tree "
                 "path:{path}\n".format(folder=tree_folder, type=starting_tree_type, path=starting_tree_path))
    curr_msa_stats["current_starting_tree_folder"] = tree_folder
    curr_msa_stats["current_starting_tree_type"] = starting_tree_type
    curr_msa_stats["starting_tree_path"] = starting_tree_path
    curr_msa_stats["tree_ind"] = i
    logging.info("Computing SPR result on full data and updating msa stats:")
    logging.info("Running naive SPR on current starting tree and updating")
    full_data_path = curr_msa_stats["local_alignment_path"]
    full_data_unique_name = curr_msa_stats["full_data_unique_name"]
    naive_spr_result_on_curr_starting_tree = SPR_analysis(full_data_path, full_data_unique_name,
                                                          curr_msa_stats,
                                                          curr_run_directory=os.path.join(curr_msa_stats.get(
                                                              "current_starting_tree_folder"), "spr_full_data_results"),
                                                          samp_lambda_function=None, full_run=True)
    logging.info(
        "SPR result on full data is {spr_result}\n Updating msa_stats".format(
            spr_result=naive_spr_result_on_curr_starting_tree))
    curr_msa_stats.update(naive_spr_result_on_curr_starting_tree)
    unique_name = "tree_{i}_lasso_spr".format(i=i)
    curr_run_directory = os.path.join(curr_msa_stats.get("current_starting_tree_folder"), unique_name)
    logging.info("Computing SPR result using lasso in folder :{}".format(curr_run_directory))
    lasso_spr_result_on_curr_starting_tree = SPR_analysis(curr_msa_stats.get("local_alignment_path"),
                                                          unique_name, curr_msa_stats, curr_run_directory,
                                                          samp_lambda_function=curr_msa_stats["lasso_predict_func"],
                                                          full_run=False)
    curr_msa_stats.update(lasso_spr_result_on_curr_starting_tree)
    complete_random_tree_stats = {k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                  k not in (
                                      IGNORE_COLS_IN_CSV
                                  )}
    logging.info("All results for starting tree {i}: \n {results}".format(i=i, results=complete_random_tree_stats))
    return complete_random_tree_stats



def rel_rf_dist_naive_vs_best_first_phase(row, curr_msa_stats, curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_first_phase"],
                                          row["naive_SPR_tree_newick"], curr_msa_stats,
                                          "naive_vs_best_first_phase",curr_run_directory)


def rel_rf_dist_naive_vs_best_second_phase(row, curr_msa_stats,curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_second_phase"],
                                          row["naive_SPR_tree_newick"], curr_msa_stats,
                                          "naive_vs_best_second_phase",curr_run_directory)


def rel_rf_dist_first_phase_vs_best(row, curr_msa_stats,curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_first_phase"],
                                          row["lasso_SPR_first_phase_tree_newick"], curr_msa_stats,
                                          "lasso_first_phase_vs_best",curr_run_directory)


def rel_rf_dist_second_phase_vs_best(row, curr_msa_stats,curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_second_phase"],
                                          row["lasso_SPR_second_phase_tree_newick"], curr_msa_stats,
                                          "lasso_second_phase_vs_best",curr_run_directory)


def get_best_ll_and_topology(ll_col, tree_topology_col):
    best_ll, max_ind = max(ll_col), ll_col.idxmax()
    best_spr_tree_newick = tree_topology_col[max_ind]
    return best_ll, best_spr_tree_newick


def enrich_curr_msa_results(curr_msa_results, curr_msa_stats,curr_run_directory):
    best_naive_spr_ll, best_naive_spr_tree_newick = get_best_ll_and_topology(
        curr_msa_results["naive_SPR_ll"],
        curr_msa_results["naive_SPR_tree_newick"])
    best_lasso_spr_first_phase_ll, best_lasso_spr_first_phase_tree_newick = get_best_ll_and_topology(
        curr_msa_results[
            "lasso_SPR_first_phase_ll"],
        curr_msa_results[
            "lasso_SPR_first_phase_tree_newick"])
    best_lasso_spr_second_phase_ll, best_lasso_spr_second_phase_tree_newick = get_best_ll_and_topology(
        curr_msa_results[
            "lasso_SPR_second_phase_ll"],
        curr_msa_results[
            "lasso_SPR_second_phase_tree_newick"])

    if best_naive_spr_ll > best_lasso_spr_first_phase_ll:
        curr_msa_results["overall_best_topology_first_phase"] = best_naive_spr_tree_newick
    else:
        curr_msa_results["overall_best_topology_first_phase"] = best_lasso_spr_first_phase_tree_newick
    if best_naive_spr_ll > best_lasso_spr_second_phase_ll:
        curr_msa_results["overall_best_topology_second_phase"] = best_naive_spr_tree_newick
    else:
        curr_msa_results["overall_best_topology_second_phase"] = best_lasso_spr_second_phase_tree_newick

    curr_msa_results["rf_naive_vs_overall_best_first_phase"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_naive_vs_best_first_phase(row, curr_msa_stats,curr_run_directory),
        axis=1)
    curr_msa_results["rf_naive_vs_overall_best_second_phase"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_naive_vs_best_second_phase(row, curr_msa_stats,curr_run_directory),
        axis=1)
    curr_msa_results["rf_first_phase_vs_overall_best"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_first_phase_vs_best(row, curr_msa_stats,curr_run_directory),
        axis=1)
    curr_msa_results["rf_second_phase_vs_overall_best"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_second_phase_vs_best(row, curr_msa_stats,curr_run_directory),
        axis=1)
    return curr_msa_results


def calculate_relative_rf_distance(best_topology_newick, given_topology_newick, curr_msa_stats,name, curr_run_directory):
    curr_run_directory = os.path.join(curr_run_directory, name)
    create_or_clean_dir(curr_run_directory)
    rf_path = os.path.join(curr_run_directory, "rf_trees_file")
    with open(rf_path, 'a+') as f_combined:
        f_combined.write(best_topology_newick)
        f_combined.write(given_topology_newick)
    relative_rf_dist = calculate_rf_dist(rf_path, curr_run_directory)
    return relative_rf_dist


def add_curr_MSA_results(n_random_starting_trees, curr_msa_stats, curr_job_output_csv_path, all_MSA_results,curr_run_directory):
    curr_msa_results = pd.DataFrame(
    )
    for i in range(n_random_starting_trees):
        random_tree_folder = os.path.join(curr_run_directory, "RANDOM_starting_tree_" + str(i))
        random_tree_path_prefix = os.path.join(random_tree_folder, "starting_tree")
        create_or_clean_dir(random_tree_folder)
        starting_tree_path = generate_random_tree_topology(curr_msa_stats["alpha"], curr_msa_stats["local_alignment_path"],
                                                           random_tree_path_prefix)
        curr_random_tree_raw_results = basic_pipeline_for_curr_starting_tree(curr_msa_stats, i, "RANDOM",
                                                                             random_tree_folder,
                                                                             starting_tree_path)
        curr_msa_results = curr_msa_results.append(curr_random_tree_raw_results, ignore_index=True)
    curr_msa_results = enrich_curr_msa_results(curr_msa_results, curr_msa_stats,curr_run_directory)
    all_MSA_results = all_MSA_results.append(curr_msa_results)
    all_MSA_results.to_csv(curr_job_output_csv_path, index=False)  # Updating data 
    return all_MSA_results



def extract_and_update_RaxML_statistics_from_full_data(curr_msa_stats):
    logging.info("Running RaxML statistics from full data and extracting statistics")
    full_data_path = curr_msa_stats["local_alignment_path"]
    full_data_unique_name = curr_msa_stats["full_data_unique_name"]
    curr_run_directory = os.path.join(curr_msa_stats.get("curr_msa_version_folder"), "raxml_full_data_results_" + \
                                      curr_msa_stats["file_name"])
    if os.path.exists(curr_run_directory):
        delete_dir_content(curr_run_directory)
    else:
        os.mkdir(curr_run_directory)
    extract_raxml_statistics_from_msa(full_data_path, full_data_unique_name, curr_msa_stats, curr_run_directory)


def get_positions_stats(alignment_df, n_seq):
    counts_per_position = [dict(alignment_df[col].value_counts()) for col in list(alignment_df)]
    gap_positions_pct = np.mean([counts_per_position[col].get('-', 0)/n_seq for col in range(len(counts_per_position))])
    position_mode = [
        sorted(counts_per_position[col].keys(), key=lambda val: counts_per_position[col][val], reverse=True)[0] for col
        in range(len(counts_per_position))]
    mode_freq_per_position = [counts_per_position[col][mode] / n_seq for col, mode in enumerate(position_mode)]
    informative_columns_count = len([i for i in range(len(counts_per_position)) if
                                     mode_freq_per_position[i] < (1 - (1 / (n_seq)))])
    probabilities = [list(map(lambda x: x / n_seq, counts_per_position[col].values())) for col in
                     list(alignment_df)]
    entropy = [sum(list(map(lambda x: -x * np.log(x), probabilities[col]))) for col in list(alignment_df)]
    avg_entropy = np.mean(entropy)
    return informative_columns_count, avg_entropy,gap_positions_pct


def generate_msa_general_stats(original_alignment_path, file_ind, current_job_results_folder, job, max_n_seq,
                               n_random_starting_trees,random_trees_training_size, random_trees_test_size
                               ):
    dataset_id = original_alignment_path
    file_name = str(file_ind)
    curr_msa_version_folder = os.path.join(current_job_results_folder, str(file_name))
    create_or_clean_dir(curr_msa_version_folder)
    logging.info("file name is " + file_name)
    full_data_unique_name = file_name
    file_type = extract_file_type(original_alignment_path, False)
    file_type_biopython = extract_file_type(original_alignment_path, True)
    with open(original_alignment_path) as original:
        original_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    orig_n_seq = len(original_alignment_data)
    local_full_msa_path = os.path.join(curr_msa_version_folder, file_name + file_type)
    take_up_to_x_sequences(original_alignment_data, local_full_msa_path, max_n_seq, file_type_biopython)
    with open(local_full_msa_path) as original:
       reduced_local_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    reduced_local_alignment_df = alignment_list_to_df(reduced_local_alignment_data )
    n_seq, n_loci = reduced_local_alignment_df.shape
    informative_columns_count, avg_entropy,gap_positions_pct = get_positions_stats(reduced_local_alignment_df, n_seq)
    curr_msa_stats = {"job_id": job, "n_seq": n_seq, "MSA_original_n_seq": orig_n_seq, "n_seq_before_reduction_by_RaxML": n_seq,
                      "n_loci": n_loci, "MSA_original_n_loci": n_loci, "dataset_id": dataset_id,
                      "curr_msa_version_folder": curr_msa_version_folder,
                      "curr_job_results_folder": current_job_results_folder,
                      "alignment_data": reduced_local_alignment_data,
                      "MSA_original_alignment_data": original_alignment_data,
                      "local_alignment_path": local_full_msa_path,
                      "full_data_unique_name": full_data_unique_name,
                      "file_ind": file_ind,
                      "file_name": file_name,
                      "file_type": file_type,
                      "file_type_biopython": file_type_biopython,
                      "n_random_starting_trees": n_random_starting_trees,
                      "random_trees_training_size": random_trees_training_size,
                      "random_trees_test_size": random_trees_test_size,
                      "max_number_of_msa_sequences": max_n_seq,
                      "informative_columns_count": informative_columns_count,
                      "avg_entropy": avg_entropy,
                      "gap_pct" : gap_positions_pct

                      }
    logging.info("Basic MSA stats computed:\n {curr_msa_stats}".format(
        curr_msa_stats={key: curr_msa_stats[key] for key in curr_msa_stats if key not in IGNORE_COLS_IN_CSV}))
    return curr_msa_stats


def re_run_on_reduced_version(curr_msa_stats, original_alignment_path, file_ind):
    raxml_reduced_file = curr_msa_stats["orig_reduced_file_path"]
    reduced_dir, rediced_fname = os.path.split(raxml_reduced_file)
    raxml_reduced_file_renamed = os.path.join(reduced_dir, curr_msa_stats["file_name"] + "_fixed" + extract_file_type(
        original_alignment_path))
    file_name_reduced = str(file_ind) + "_fixed"
    os.rename(raxml_reduced_file, raxml_reduced_file_renamed)
    logging.warning("Reduced version of previous file is found in " + raxml_reduced_file_renamed)
    logging.warning("Re calculating MSA stats on reduced version")
    with open(raxml_reduced_file_renamed) as reduced_path:
        file_type_biopython = extract_file_type(raxml_reduced_file_renamed, True)
        reduced_data = list(SeqIO.parse(reduced_path, file_type_biopython))
        n_loci_reduced = len(reduced_data[0].seq)
        n_seq_reduced = len(reduced_data)  # supposed to stay the same as "n_seq_before_reduction"
    original_alignment_df_reduced = alignment_list_to_df(reduced_data)
    informative_columns_count_reduced, avg_entropy_reduced,gap_pct_reduced = get_positions_stats(original_alignment_df_reduced,
                                                                                 n_seq_reduced)
    reduced_curr_msa_stats = curr_msa_stats.copy()
    reduced_curr_msa_stats.update(
        {"file_name": file_name_reduced, "local_alignment_path": raxml_reduced_file_renamed, "n_seq": n_seq_reduced,
         "n_loci": n_loci_reduced,
         "alignment_data": reduced_data,
         "informative_columns_count": informative_columns_count_reduced,
         "avg_entropy": avg_entropy_reduced,
         "gap_pct" : gap_pct_reduced
         })
    logging.info("New msa stats for reduced data: {msa_stats}".format(msa_stats=reduced_curr_msa_stats))
    return reduced_curr_msa_stats




def main():
    args = job_parser()
    job_ind, curr_job_folder, max_n_sequences, n_random_starting_trees, random_trees_training_size, random_trees_test_size, only_evaluate_lasso = args.job_ind, args.curr_job_folder, args.max_n_sequences, \
                                                                                                                          args.n_random_starting_trees, args.random_trees_training_size,args.random_trees_test_size, args.only_evaluate_lasso
    job_related_file_paths = get_job_related_files_paths(curr_job_folder, job_ind)
    job_msa_paths_file, general_log_path, job_csv_path, spr_log_path, curr_job_status_file = job_related_file_paths[
                                                                                                 "job_msa_paths_file"], \
                                                                                             job_related_file_paths[
                                                                                                 "general_log_path"], \
                                                                                             job_related_file_paths[
                                                                                                 "job_csv_path"], \
                                                                                             job_related_file_paths[
                                                                                                 "spr_log_path"], \
                                                                                             job_related_file_paths[
                                                                                                 "job_status_file"]
    with open(job_msa_paths_file, "r") as paths_file:
        curr_job_file_path_list = paths_file.read().splitlines()
    logging.basicConfig(filename=general_log_path, level=LOGGING_LEVEL)
    logging.info('#Started running on job' + str(job_ind))
    all_MSA_results = pd.DataFrame(
    )
    all_MSA_results.to_csv(job_csv_path, index=False)
    for file_ind, original_alignment_path in enumerate(curr_job_file_path_list):
        logging.info(' #running on file ind ' + str(file_ind) + " path=" + str(original_alignment_path))
        curr_msa_stats = generate_msa_general_stats(
            original_alignment_path, file_ind, curr_job_folder, job_ind, max_n_sequences, n_random_starting_trees,random_trees_training_size, random_trees_test_size)
        try:

            logging.info("Computing raxml result on full data:")
            extract_and_update_RaxML_statistics_from_full_data(
                curr_msa_stats)
        except RE_RUN_ON_REDUCED_VERSION:
            curr_msa_stats = re_run_on_reduced_version(curr_msa_stats, original_alignment_path,
                                                       file_ind
                                                       )  # use current curr_msa_stats
            extract_and_update_RaxML_statistics_from_full_data(
                curr_msa_stats)
        curr_msa_stats["spr_log_path"] = spr_log_path
        curr_msa_stats["current_job_results_folder"] = curr_job_folder
        brlen_generators = {'exponential':sample_exp,'uniform': sample_uniform,'optimized': None}
        if random_trees_training_size==-1:
            training_size_options = TRAINING_SIZE_OPTIONS
        else:
            training_size_options = [random_trees_training_size]
        for brlen_generator_name in  brlen_generators:
            brlen_run_directory = os.path.join(curr_msa_stats["curr_msa_version_folder"], brlen_generator_name)
            create_dir_if_not_exists(brlen_run_directory)
            for training_size in training_size_options:
                curr_run_directory = os.path.join(brlen_run_directory, str(training_size))
                create_dir_if_not_exists(curr_run_directory)
                comb_name = "{brlen_generator_name}_s_{training_size}".format(
                    brlen_generator_name=brlen_generator_name, training_size=training_size)
                logging.info("Using {} to generate branch lengths".format(brlen_generator_name))
                brlen_generator_func =brlen_generators.get(brlen_generator_name)
                curr_msa_stats["brlen_generator"] = brlen_generator_name
                curr_msa_stats["actucal_training_size"] = training_size
                logging.info("Generation training dataset in folder: {}".format(curr_run_directory))
                training_sitelh = generate_site_lh_data(curr_msa_stats=curr_msa_stats,n_iter=training_size,name="training_"+ comb_name,
                                                        brlen_generator_func=brlen_generator_func, curr_run_directory=curr_run_directory)
                curr_msa_stats["training_sitelh_df"] = training_sitelh
                logging.info("Generation test dataset in folder: {}".format(curr_run_directory))
                test_sitelh = generate_site_lh_data(curr_msa_stats=curr_msa_stats, n_iter=random_trees_test_size,
                                                        name="test_"+ comb_name,
                                                        brlen_generator_func=brlen_generator_func,curr_run_directory=curr_run_directory)
                curr_msa_stats["test_sitelh_df"] = test_sitelh
                apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats,name=brlen_generator_name, curr_run_directory=curr_run_directory)  # calculating positions_weight
                if only_evaluate_lasso:
                    logging.info("only evaluating lasso: " + str(curr_msa_stats))
                    lasso_evaluation_result = {k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                               k not in IGNORE_COLS_IN_CSV
                                               }
                    all_MSA_results = all_MSA_results.append(lasso_evaluation_result, ignore_index=True)
                    all_MSA_results.to_csv(job_csv_path, index=False)
                else:
                    all_MSA_results = add_curr_MSA_results(n_random_starting_trees, curr_msa_stats, job_csv_path,
                                                           all_MSA_results,curr_run_directory)
    with open(curr_job_status_file, 'w') as job_status_f:
        job_status_f.write("Done")



if __name__ == "__main__":
    main()
