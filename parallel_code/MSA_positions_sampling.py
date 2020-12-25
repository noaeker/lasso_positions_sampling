from lasso_model_pipeline import *
from generate_SPR import *
from sys import argv
from shutil import copyfile
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


def find_best_topology_up_to_first_phase(row):
    if row["best_naive_spr_ll"] >= row["best_lasso_SPR_first_phase_ll"]:
        val = row["best_naive_spr_tree_topology_newick"]
    else:
        val = row["best_lasso_SPR_first_phase_tree_newick"]
    return val


def find_ll_winner_up_to_first_phase(row):
    if row["best_naive_spr_ll"] >= row["best_lasso_SPR_first_phase_ll"]:
        val = 'naive_spr'
    else:
        val = 'lasso_first_spr'
    return val


def find_best_topology_up_to_second_phase(row):
    if row["best_naive_spr_ll"] >= row["best_lasso_SPR_second_phase_ll"]:
        val = row["best_naive_spr_tree_topology_newick"]
    else:
        val = row["best_lasso_SPR_second__phase_tree_newick"]
    return val


def find_ll_winner_up_to_second_phase(row):
    if row["best_naive_spr_ll"] >= row["best_lasso_SPR_second_phase_ll"]:
        val = 'naive_spr'
    else:
        val = 'lasso_second_spr'
    return val


def calc_relative_rf_distance_up_to_first_phase(row, curr_msa_stats):
    return calculate_relative_rf_distance(row["best_naive_spr_tree_topology_newick"],
                                          row["best_lasso_SPR_first_phase_tree_newick"], curr_msa_stats,
                                          "naive_vs_first_phase")


def calc_relative_rf_distance_up_to_second_phase(row, curr_msa_stats):
    return calculate_relative_rf_distance(row["best_naive_spr_tree_topology_newick"],
                                          row["best_lasso_SPR_second_phase_tree_newick"], curr_msa_stats,
                                          "naive_vs_second_phase")


def enrich_curr_msa_results(curr_msa_results, curr_msa_stats):
    curr_msa_results["best_naive_spr_ll"], max_index_naive_spr = max(curr_msa_results["naive_SPR_ll"]), \
                                                                 curr_msa_results["naive_SPR_ll"].idxmax()
    curr_msa_results["best_naive_spr_tree_topology_newick"] = curr_msa_results["naive_SPR_tree_newick"][
        max_index_naive_spr]

    curr_msa_results["best_lasso_SPR_first_phase_ll"], max_index_lasso_first_phase = max(
        curr_msa_results["lasso_SPR_first_phase_ll"]), \
                                                                                     curr_msa_results[
                                                                                         "lasso_SPR_first_phase_ll"].idxmax()
    curr_msa_results["best_lasso_SPR_first_phase_tree_newick"] = curr_msa_results["lasso_SPR_first_phase_tree_newick"][
        max_index_lasso_first_phase]

    curr_msa_results["best_lasso_SPR_second_phase_ll"], max_index_lasso_second_phase = max(
        curr_msa_results["lasso_SPR_second_phase_ll"]), \
                                                                                       curr_msa_results[
                                                                                           "lasso_SPR_first_phase_ll"].idxmax()
    curr_msa_results["best_lasso_SPR_second_phase_tree_newick"] = \
        curr_msa_results["lasso_SPR_second_phase_tree_newick"][
            max_index_lasso_second_phase]
    curr_msa_results["ll_winner_up_to_first_phase"] = curr_msa_results.apply(find_ll_winner_up_to_first_phase, axis=1)
    curr_msa_results["best_topology_up_to_first_phase"] = curr_msa_results.apply(find_best_topology_up_to_first_phase,
                                                                                 axis=1)
    curr_msa_results["ll_winner_up_to_second_phase"] = curr_msa_results.apply(find_ll_winner_up_to_second_phase, axis=1)
    curr_msa_results["best_topology_up_to_second_phase"] = curr_msa_results.apply(find_best_topology_up_to_second_phase,
                                                                                  axis=1)
    curr_msa_results["rf_naive_vs_first_phase"] = curr_msa_results.apply(
        lambda row: calc_relative_rf_distance_up_to_first_phase(row, curr_msa_stats),
        axis=1)
    curr_msa_results["rf_naive_vs_second_phase"] = curr_msa_results.apply(
        lambda row: calc_relative_rf_distance_up_to_second_phase(row, curr_msa_stats),
        axis=1)
    return curr_msa_results


def calculate_relative_rf_distance(best_topology_newick, given_topology_newick, curr_msa_stats, name):
    curr_run_directory = os.path.join(curr_msa_stats.get(
        "curr_msa_version_folder"), name)
    create_or_clean_dir(curr_run_directory)
    rf_path = os.path.join(curr_run_directory, "rf_trees_file")
    with open(rf_path, 'a+') as f_combined:
        f_combined.write(best_topology_newick)
        f_combined.write(given_topology_newick)
    relative_rf_dist = calculate_rf_dist(rf_path, curr_run_directory)
    return relative_rf_dist


def add_curr_MSA_results(n_random_starting_trees, curr_msa_stats, curr_job_output_csv_path, all_MSA_results):
    curr_msa_results = pd.DataFrame(
    )
    for i in range(n_random_starting_trees):
        random_tree_folder = os.path.join(curr_msa_stats.get(
            "curr_msa_version_folder"), "RANDOM_starting_tree_" + str(i))
        random_tree_path_prefix = os.path.join(random_tree_folder, "starting_tree")
        create_or_clean_dir(random_tree_folder)
        starting_tree_path = generate_random_tree(curr_msa_stats["alpha"], curr_msa_stats["local_alignment_path"],
                                                  random_tree_path_prefix)
        curr_random_tree_raw_results = basic_pipeline_for_curr_starting_tree(curr_msa_stats, i, "RANDOM",
                                                                             random_tree_folder,
                                                                             starting_tree_path)
        curr_msa_results = curr_msa_results.append(curr_random_tree_raw_results, ignore_index=True)
    curr_msa_results = enrich_curr_msa_results(curr_msa_results, curr_msa_stats)
    all_MSA_results = all_MSA_results.append(curr_msa_results)
    all_MSA_results.to_csv(curr_job_output_csv_path, index=False)  # Updating data after each tree
    return all_MSA_results


def generate_site_lh_data(curr_msa_stats, n_iter):
    curr_sitelh_folder = os.path.join(curr_msa_stats["curr_msa_version_folder"], "sitelh")
    create_dir_if_not_exists(curr_sitelh_folder)
    curr_csv_path = os.path.join(curr_sitelh_folder, curr_msa_stats.get("file_name") + ".csv")
    local_file_path = curr_msa_stats.get("local_alignment_path")
    logging.info(
        "Computing sitelh data on " + local_file_path + " from beggining,since no such csv path found " + curr_csv_path)
    logging.info("Generating " + str(n_iter) + " random trees ")
    sitelh_ll_list = []
    for i in range(n_iter):
        curr_site_lh = raxml_compute_per_site_ll_on_a_random_tree(local_file_path, i, curr_msa_stats)
        sitelh_ll_list.append(curr_site_lh)
    sitelh_df = pd.DataFrame(sitelh_ll_list, columns=list(range(len(sitelh_ll_list[0]))),
                             index=list(range(len(sitelh_ll_list))))
    logging.info(
        "Writing sitelh data to  " + curr_csv_path)
    sitelh_df.to_csv(
        curr_csv_path, index=False)
    logging.info(
        "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=curr_csv_path))
    curr_msa_stats["sitelh_df"] = sitelh_df
    return sitelh_df


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
    run_raxml_on_full_dataset(full_data_path, full_data_unique_name, curr_msa_stats, curr_run_directory)


def get_positions_stats(original_alignment_df, n_seq):
    counts_per_position = [dict(original_alignment_df[col].value_counts()) for col in list(original_alignment_df)]
    position_mode = [
        sorted(counts_per_position[col].keys(), key=lambda val: counts_per_position[col][val], reverse=True)[0] for col
        in range(len(counts_per_position))]
    mode_freq_per_position = [counts_per_position[col][mode] / n_seq for col, mode in enumerate(position_mode)]
    informative_columns_count = len([i for i in range(len(counts_per_position)) if
                                     mode_freq_per_position[i] < (1 - (1 / (n_seq)))])
    probabilities = [list(map(lambda x: x / n_seq, counts_per_position[col].values())) for col in
                     list(original_alignment_df)]
    entropy = [sum(list(map(lambda x: -x * np.log(x), probabilities[col]))) for col in list(original_alignment_df)]
    avg_entropy = np.median(entropy)
    return informative_columns_count, avg_entropy


def generate_msa_general_stats(original_alignment_path, file_ind, current_job_results_folder, job, max_n_seq,
                               n_random_starting_trees
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
        original_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    original_alignment_df = alignment_list_to_df(original_alignment_data)
    informative_columns_count, avg_entropy = get_positions_stats(original_alignment_df, orig_n_seq)
    n_seq, n_loci = original_alignment_df.shape
    curr_msa_stats = {"job_id": job, "n_seq": n_seq, "original_n_seq": orig_n_seq, "n_seq_before_reduction": n_seq,
                      "n_loci": n_loci, "original_n_loci": n_loci, "dataset_id": dataset_id,
                      "curr_msa_version_folder": curr_msa_version_folder,
                      "curr_job_results_folder": current_job_results_folder,
                      "alignment_data": original_alignment_data,
                      "local_alignment_path": local_full_msa_path,
                      "full_data_unique_name": full_data_unique_name,
                      "file_ind": file_ind,
                      "file_name": file_name,
                      "file_type": file_type,
                      "file_type_biopython": file_type_biopython,
                      "random_trees_sample_size": n_random_starting_trees,
                      "max_number_of_msa_sequences": max_n_seq,
                      "informative columns count": informative_columns_count,
                      "avg_entropy": avg_entropy

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
    informative_columns_count_reduced, avg_entropy_reduced = get_positions_stats(original_alignment_df_reduced,
                                                                                 n_seq_reduced)
    reduced_curr_msa_stats = curr_msa_stats.copy()
    reduced_curr_msa_stats.update(
        {"file_name": file_name_reduced, "local_alignment_path": raxml_reduced_file_renamed, "n_seq": n_seq_reduced,
         "n_loci": n_loci_reduced,
         "alignment_data": reduced_data,
         "informative columns count": informative_columns_count_reduced,
         "avg_entropy": avg_entropy_reduced
         })
    logging.info("New msa stats for reduced data: {msa_stats}".format(msa_stats=reduced_curr_msa_stats))
    return reduced_curr_msa_stats


def main():
    args = job_parser()
    job_ind, curr_job_folder, max_n_sequences, n_random_starting_trees, only_evaluate_lasso = args.job_ind, args.curr_job_folder, args.max_n_sequences, \
                                                                                              args.n_random_starting_trees, args.only_evaluate_lasso
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
    for file_ind, original_alignment_path in enumerate(curr_job_file_path_list):
        all_MSA_results = pd.DataFrame(
        )
        all_MSA_results.to_csv(job_csv_path, index=False)
        logging.info(' #running on file ind ' + str(file_ind) + " path=" + str(original_alignment_path))
        curr_msa_stats = generate_msa_general_stats(
            original_alignment_path, file_ind, curr_job_folder, job_ind, max_n_sequences, n_random_starting_trees)
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
        generate_site_lh_data(curr_msa_stats, RANDOM_TREES_TRAINING_SIZE)
        # curr_msa_stats["raxml_parsimony_tree_path"]
        apply_lasso_on_data_and_update_stats(curr_msa_stats)  # calculating positions_weight
        if only_evaluate_lasso:
            logging.info("only evaluating lasso: " + str(curr_msa_stats))
            lasso_evaluation_result = {k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                       k not in IGNORE_COLS_IN_CSV
                                       }
            all_MSA_results = all_MSA_results.append(lasso_evaluation_result, ignore_index=True)
            all_MSA_results.to_csv(job_csv_path, index=False)
        else:
            all_MSA_results = add_curr_MSA_results(n_random_starting_trees, curr_msa_stats, job_csv_path,
                                                   all_MSA_results)
        with open(curr_job_status_file, 'w') as job_status_f:
            job_status_f.write("Done")
            # delete_file_content(general_log_path)
            # delete_file_content(spr_log_path)


if __name__ == "__main__":
    main()
