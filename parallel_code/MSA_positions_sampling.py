from lasso_model_pipeline import *
from generate_SPR import *
from training_and_test_set_generation import *
from raxml import *
import filecmp


def rel_rf_dist_naive_vs_best_first_phase(row, curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_first_phase"],
                                          row["naive_SPR_tree_newick"],
                                          "naive_vs_best_first_phase", curr_run_directory)


def rel_rf_dist_naive_vs_best_second_phase(row, curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_second_phase"],
                                          row["naive_SPR_tree_newick"],
                                          "naive_vs_best_second_phase", curr_run_directory)


def rel_rf_dist_first_phase_vs_best(row, curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_first_phase"],
                                          row["lasso_SPR_first_phase_tree_newick"],
                                          "lasso_first_phase_vs_best", curr_run_directory)


def rel_rf_dist_second_phase_vs_best(row, curr_run_directory):
    return calculate_relative_rf_distance(row["overall_best_topology_second_phase"],
                                          row["lasso_SPR_second_phase_tree_newick"],
                                          "lasso_second_phase_vs_best", curr_run_directory)


def get_best_ll_and_topology(ll_col, tree_topology_col):
    best_ll, max_ind = max(ll_col), ll_col.idxmax()
    best_spr_tree_newick = tree_topology_col[max_ind]
    return best_ll, best_spr_tree_newick


def enrich_curr_msa_results(curr_msa_results, curr_run_directory):
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
        lambda row: rel_rf_dist_naive_vs_best_first_phase(row, curr_run_directory),
        axis=1)
    curr_msa_results["rf_naive_vs_overall_best_second_phase"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_naive_vs_best_second_phase(row, curr_run_directory),
        axis=1)
    curr_msa_results["rf_first_phase_vs_overall_best"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_first_phase_vs_best(row, curr_run_directory),
        axis=1)
    curr_msa_results["rf_second_phase_vs_overall_best"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_second_phase_vs_best(row, curr_run_directory),
        axis=1)
    return curr_msa_results


def calculate_relative_rf_distance(best_topology_newick, given_topology_newick, name,
                                   curr_run_directory):
    curr_run_directory = os.path.join(curr_run_directory, name)
    create_or_clean_dir(curr_run_directory)
    rf_path = os.path.join(curr_run_directory, "rf_trees_file")
    with open(rf_path, 'a+') as f_combined:
        f_combined.write(best_topology_newick)
        f_combined.write(given_topology_newick)
    relative_rf_dist = calculate_rf_dist(rf_path, curr_run_directory)
    return relative_rf_dist


def generate_random_starting_tree(i, curr_run_directory, curr_msa_stats):
    random_tree_folder = os.path.join(curr_run_directory, "RANDOM_starting_tree_" + str(i))
    random_tree_path_prefix = os.path.join(random_tree_folder, "starting_tree")
    create_or_clean_dir(random_tree_folder)
    starting_tree_path = random_tree_path_prefix + ".raxml.startTree"
    baseline_starting_tree_path = starting_tree_path.replace(curr_msa_stats["run_prefix"],
                                                             curr_msa_stats["baseline_run_prefix"])
    if os.path.exists(baseline_starting_tree_path):
        shutil.copyfile(baseline_starting_tree_path, starting_tree_path)
    if not os.path.exists(starting_tree_path):
        logging.info("Generating a totally random tree as a starting tree")
        starting_tree_path = generate_random_tree_topology(curr_msa_stats["alpha"],
                                                           curr_msa_stats["local_alignment_path"],
                                                           random_tree_path_prefix)
    return starting_tree_path


def naive_spr_on_current_starting_tree(starting_tree_path, curr_msa_stats, curr_run_directory):
    logging.info("Running naive SPR on current starting tree and updating")
    full_data_path = curr_msa_stats["local_alignment_path"]
    curr_run_directory = os.path.join(curr_run_directory, "spr_full_data_results")
    naive_spr_result_on_curr_starting_tree = SPR_analysis(full_data_path, starting_tree_path,
                                                          curr_msa_stats,
                                                          curr_run_directory=curr_run_directory,
                                                          full_run=True)
    return naive_spr_result_on_curr_starting_tree


def lasso_based_spr_on_current_starting_tree(starting_tree_path, curr_msa_stats, curr_run_directory):
    lasso_spr_result_on_curr_starting_tree = SPR_analysis(curr_msa_stats.get("local_alignment_path"),
                                                          starting_tree_path, curr_msa_stats,
                                                          curr_run_directory=curr_run_directory, full_run=False)
    return lasso_spr_result_on_curr_starting_tree


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
    gap_positions_pct = np.mean(
        [counts_per_position[col].get('-', 0) / n_seq for col in range(len(counts_per_position))])
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
    return informative_columns_count, avg_entropy, gap_positions_pct


def generate_msa_general_stats(original_alignment_path, file_ind, current_msa_version_folder,
                               current_job_results_folder, job, max_n_seq,
                               n_random_starting_trees, random_trees_training_size, random_trees_test_size, run_prefix,
                               baseline_run_prefix
                               ):
    dataset_id = original_alignment_path
    file_name = str(file_ind)
    logging.info("file name is " + file_name)
    full_data_unique_name = file_name
    file_type = extract_file_type(original_alignment_path, False)
    file_type_biopython = extract_file_type(original_alignment_path, True)
    with open(original_alignment_path) as original:
        original_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    orig_n_seq = len(original_alignment_data)
    local_full_msa_path = os.path.join(current_msa_version_folder, file_name + file_type)
    take_up_to_x_sequences(original_alignment_data, local_full_msa_path, max_n_seq, file_type_biopython)
    with open(local_full_msa_path) as original:
        reduced_local_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    reduced_local_alignment_df = alignment_list_to_df(reduced_local_alignment_data)
    n_seq, n_loci = reduced_local_alignment_df.shape
    informative_columns_count, avg_entropy, gap_positions_pct = get_positions_stats(reduced_local_alignment_df, n_seq)
    curr_msa_stats = {"job_id": job, "n_seq": n_seq, "MSA_original_n_seq": orig_n_seq,
                      "n_seq_before_reduction_by_RaxML": n_seq,
                      "n_loci": n_loci, "MSA_original_n_loci": n_loci, "dataset_id": dataset_id,
                      "curr_msa_version_folder": current_msa_version_folder,
                      "curr_job_results_folder": current_job_results_folder,
                      "run_prefix": run_prefix,
                      "baseline_run_prefix": baseline_run_prefix,
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
                      "gap_pct": gap_positions_pct

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
    informative_columns_count_reduced, avg_entropy_reduced, gap_pct_reduced = get_positions_stats(
        original_alignment_df_reduced,
        n_seq_reduced)
    reduced_curr_msa_stats = curr_msa_stats.copy()
    update_dict = {"file_name": file_name_reduced, "local_alignment_path": raxml_reduced_file_renamed,
                   "n_seq": n_seq_reduced,
                   "n_loci": n_loci_reduced,
                   "alignment_data": reduced_data,
                   "informative_columns_count": informative_columns_count_reduced,
                   "avg_entropy": avg_entropy_reduced,
                   "gap_pct": gap_pct_reduced
                   }
    reduced_curr_msa_stats.update(
        update_dict)
    logging.info(
        "New msa stats for reduced data: {msa_stats}".format(msa_stats=update_dict.copy().pop("alignment_data")))
    return reduced_curr_msa_stats


def update_chosen_brlen_generators(exp_brlen, uni_brlen, opt_brlen):
    brlen_generators = {}
    for brlen_generator_name, brlen_generator_func in BRLEN_GENERATORS.items():
        if (brlen_generator_name == 'exponential' and exp_brlen) or (
                brlen_generator_name == 'uniform' and uni_brlen) or (
                brlen_generator_name == 'optimized' and opt_brlen):
            brlen_generators[brlen_generator_name] = brlen_generator_func
    return brlen_generators
    logging.info("Brlen genertors are chosen to be {}".format(str(brlen_generators.keys())))




def Lasso_training_and_test(brlen_generators, curr_msa_stats, training_size_options, random_trees_test_size):
    Lasso_folder = os.path.join(curr_msa_stats["curr_msa_version_folder"], "Lasso_folder")
    create_dir_if_not_exists(Lasso_folder)
    logging.info("Generating Lasso folder in {}".format(Lasso_folder))
    random_trees_folder = os.path.join(Lasso_folder, "random_tree_generation")
    curr_msa_stats["Lasso_folder"] = Lasso_folder
    curr_msa_stats["random_trees_folder"] = random_trees_folder
    create_dir_if_not_exists(random_trees_folder)
    random_trees_per_training_size = {}
    for training_size in training_size_options:
        training_random_tree_path_and_folder_list = generate_n_random_topologies(training_size, random_trees_folder,
                                                                                 curr_msa_stats, "training")
        random_trees_per_training_size[training_size] = training_random_tree_path_and_folder_list
    random_trees_test = generate_n_random_topologies(random_trees_test_size, random_trees_folder, curr_msa_stats,
                                                     "test")
    run_configurations = {}
    for brlen_generator_name in brlen_generators:
        brlen_run_directory = os.path.join(Lasso_folder, brlen_generator_name)
        create_dir_if_not_exists(brlen_run_directory)
        brlen_generator_func = brlen_generators.get(brlen_generator_name)
        test_folder = os.path.join(brlen_run_directory, "test_{}_random_trees_eval".format(random_trees_test_size))
        create_dir_if_not_exists(test_folder)
        test_sitelh, test_sitelh_path = get_test_set_df(curr_msa_stats, brlen_generator_func,
                                                        test_folder, random_trees_test)
        for training_size in training_size_options:
            training_size_directory = os.path.join(brlen_run_directory,
                                                   "training_{}_random_tree_eval".format(training_size))
            create_dir_if_not_exists(training_size_directory)
            logging.info("")
            training_sitelh, training_sitelh_path = get_training_df(curr_msa_stats, brlen_generator_func,
                                                                    training_size_directory,
                                                                    random_trees_per_training_size[training_size])
            Lasso_results = apply_lasso_on_sitelh_data_and_update_statistics(curr_msa_stats,
                                                                             curr_run_directory=training_size_directory,
                                                                             sitelh_training_df=training_sitelh,
                                                                             sitelh_test_df=test_sitelh)  # calculating positions_weight
            if brlen_generator_name not in run_configurations:
                run_configurations[brlen_generator_name] = {}
            run_configurations[brlen_generator_name][training_size] = Lasso_results
    return run_configurations



def generate_msa_stats_and_lasso_from_beggining(original_alignment_path, file_ind, curr_msa_version_folder,
                                                curr_job_folder, job_ind, max_n_sequences,
                                                n_random_starting_trees,
                                                random_trees_training_size, random_trees_test_size, run_prefix,
                                                baseline_run_prefix, spr_log_path,
                                                brlen_generators, training_size_options):
    logging.info(
        "Generating msa stats and Lasso in {} from beggining (not using a baseline)".format(
            curr_msa_version_folder))
    curr_msa_stats = generate_msa_general_stats(
        original_alignment_path, file_ind, curr_msa_version_folder, curr_job_folder, job_ind, max_n_sequences,
        n_random_starting_trees,
        random_trees_training_size, random_trees_test_size, run_prefix, baseline_run_prefix)
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
    curr_msa_stats_dump_path = os.path.join(curr_msa_stats["curr_msa_version_folder"], "curr_msa_stats.dump")
    lasso_configurations_per_training_size = Lasso_training_and_test(brlen_generators, curr_msa_stats,
                                                                     training_size_options,
                                                                     random_trees_test_size)
    with open(curr_msa_stats_dump_path, 'wb') as handle:
        pickle.dump(curr_msa_stats, handle, protocol=pickle.HIGHEST_PROTOCOL)
    curr_msa_stats_dump_path = os.path.join(curr_msa_stats["curr_msa_version_folder"], "curr_msa_stats.dump")
    with open(curr_msa_stats_dump_path, 'wb') as handle:
        pickle.dump(curr_msa_stats, handle, protocol=pickle.HIGHEST_PROTOCOL)
    lasso_configurations_per_training_size_dump_path = os.path.join(curr_msa_stats["curr_msa_version_folder"],
                                                                    "lasso.dump")
    with open(lasso_configurations_per_training_size_dump_path, 'wb') as handle:
        pickle.dump(lasso_configurations_per_training_size, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return curr_msa_stats, lasso_configurations_per_training_size



def perform_only_lasso_pipeline(training_size_options, brlen_generators, curr_msa_stats, lasso_configurations_per_training_size,
                                job_csv_path):
    all_msa_results = pd.DataFrame(
    )
    all_msa_results.to_csv(job_csv_path, index=False)
    for brlen_generator_name in brlen_generators:
        curr_msa_stats["brlen_generator"] = brlen_generator_name
        for training_size in training_size_options:
            curr_msa_stats["actucal_training_size"] = training_size
            lasso_results = lasso_configurations_per_training_size[brlen_generator_name][training_size]
            curr_msa_stats.update(lasso_results)
            logging.info("only evaluating lasso: ")
            lasso_evaluation_result = {k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                       k not in IGNORE_COLS_IN_CSV
                                       }
            all_msa_results = all_msa_results.append(lasso_evaluation_result, ignore_index=True)
            all_msa_results.to_csv(job_csv_path, index=False)


def perform_spr_pipeline(training_size_options, brlen_generators, curr_msa_stats, lasso_configurations_per_training_size,
                                job_csv_path):
    all_msa_results = pd.DataFrame(
    )
    all_msa_results.to_csv(job_csv_path, index=False)
    for i in range(curr_msa_stats["n_random_starting_trees"]):
        starting_tree_run_directory = os.path.join(curr_msa_stats["curr_msa_version_folder"],
                                                   'starting_tree_{i}'.format(i=i))
        create_dir_if_not_exists(starting_tree_run_directory)
        logging.info(
            "Working on starting tree {i} in dicrectory {dir}".format(i=i, dir=starting_tree_run_directory))
        starting_tree_path = generate_random_starting_tree(i, starting_tree_run_directory, curr_msa_stats)
        curr_msa_stats["starting_tree_path"] = starting_tree_path
        naive_spr_results_dump = os.path.join(starting_tree_run_directory, 'naive_spr.dump')
        naive_spr_results_dump_baseline = naive_spr_results_dump.replace(curr_msa_stats["run_prefix"],curr_msa_stats["baseline_run_prefix"])
        if os.path.exists(naive_spr_results_dump_baseline):
            with open(naive_spr_results_dump_baseline, 'rb') as handle:
                naive_spr_results = pickle.load(handle)
        else:
            naive_spr_results = naive_spr_on_current_starting_tree(starting_tree_path, curr_msa_stats,
                                                                   curr_run_directory=starting_tree_run_directory)
            with open(naive_spr_results_dump, 'wb') as handle:
                pickle.dump(naive_spr_results_dump, handle, protocol=pickle.HIGHEST_PROTOCOL)
        curr_msa_stats.update(naive_spr_results)
        for brlen_generator_name in brlen_generators:
            brlen_run_directory = os.path.join(starting_tree_run_directory, brlen_generator_name)
            create_dir_if_not_exists(brlen_run_directory)
            curr_msa_stats["brlen_generator"] = brlen_generator_name
            for training_size in training_size_options:
                curr_msa_stats["actucal_training_size"] = training_size
                curr_training_size_directory = os.path.join(brlen_run_directory, str(training_size))
                create_dir_if_not_exists(curr_training_size_directory)
                lasso_results = lasso_configurations_per_training_size[brlen_generator_name][training_size]
                curr_msa_stats.update(lasso_results)
                logging.info("Starting Lasso-based SPR on training size: {size}")
                lasso_based_spr_results = lasso_based_spr_on_current_starting_tree(starting_tree_path,
                                                                                   curr_msa_stats,
                                                                                   curr_training_size_directory)
                curr_msa_stats.update(lasso_based_spr_results)

                all_msa_results = all_msa_results.append({k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                                          k not in IGNORE_COLS_IN_CSV
                                                          }, ignore_index=True)
                all_msa_results.to_csv(job_csv_path)






def main():
    args = job_parser()
    job_ind, curr_job_folder, max_n_sequences, n_random_starting_trees, random_trees_training_size, random_trees_test_size, run_prefix, baseline_run_prefix, only_evaluate_lasso, exp_brlen, uni_brlen, opt_brlen = args.job_ind, args.curr_job_folder, args.max_n_sequences, \
                                                                                                                                                                                                                    args.n_random_starting_trees, args.random_trees_training_size, args.random_trees_test_size, args.run_prefix, args.baseline_run_prefix, args.only_evaluate_lasso, args.exp_brlen, args.uni_brlen, args.opt_brlen
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
    baseline_msa_paths_file = job_msa_paths_file.replace(run_prefix, baseline_run_prefix)
    with open(job_msa_paths_file, "r") as paths_file:
        curr_job_file_path_list = paths_file.read().splitlines()
    logging.basicConfig(filename=general_log_path, level=LOGGING_LEVEL)
    logging.info('#Started running on job' + str(job_ind))
    if not os.path.exists(baseline_msa_paths_file):
        logging.info("Not using baseline for this job")
    elif filecmp.cmp(baseline_msa_paths_file, job_msa_paths_file):
        logging.info("Files in the baseline folder are matching")
    else:
        logging.error("Files in the baseline aren't matching!!!")
    for file_ind, original_alignment_path in enumerate(curr_job_file_path_list):
        msa_name = original_alignment_path.replace(MSAs_FOLDER, "").replace("ref_msa.aa.phy", "").replace(os.path.sep,
                                                                                                          "_")
        logging.info('#running on file name {} and ind (relativ to job) {} path= {}'.format(msa_name, file_ind,
                                                                                            original_alignment_path))
        curr_msa_version_folder = os.path.join(curr_job_folder, msa_name)
        create_or_clean_dir(curr_msa_version_folder)
        curr_msa_version_lasso_dump = os.path.join(curr_msa_version_folder, 'lasso.dump')
        curr_msa_version_stats_dump = os.path.join(curr_msa_version_folder, 'curr_msa_stats.dump')
        curr_msa_version_folder_baseline = curr_msa_version_folder.replace(run_prefix, baseline_run_prefix)
        curr_msa_version_lasso_dump_baseline = curr_msa_version_lasso_dump.replace(run_prefix, baseline_run_prefix)
        curr_msa_version_stats_dump_baseline = curr_msa_version_stats_dump.replace(run_prefix, baseline_run_prefix)
        brlen_generators = update_chosen_brlen_generators(exp_brlen, uni_brlen, opt_brlen)
        if random_trees_training_size == -1:
            training_size_options = TRAINING_SIZE_OPTIONS
        else:
            training_size_options = [random_trees_training_size]
        if os.path.exists(curr_msa_version_folder_baseline) and os.path.exists(
                curr_msa_version_lasso_dump_baseline) and os.path.exists(curr_msa_version_stats_dump_baseline):
            logging.info(
                "Using dump files in {} to generate msa stats and Lasso".format(curr_msa_version_folder_baseline))
            with open(curr_msa_version_stats_dump_baseline, 'rb') as handle:
                curr_msa_stats = pickle.load(handle)
                curr_msa_stats["curr_msa_version_folder"] = curr_msa_version_folder
            with open(curr_msa_version_lasso_dump_baseline, 'rb') as handle:
                lasso_configurations_per_training_size = pickle.load(handle)

        else:
            curr_msa_stats, lasso_configurations_per_training_size = generate_msa_stats_and_lasso_from_beggining(
                original_alignment_path, file_ind, curr_msa_version_folder, curr_job_folder, job_ind, max_n_sequences,
                n_random_starting_trees,
                random_trees_training_size, random_trees_test_size, run_prefix, baseline_run_prefix, spr_log_path,
                brlen_generators, training_size_options)
        if only_evaluate_lasso:
            perform_only_lasso_pipeline(training_size_options, brlen_generators, curr_msa_stats,
                                        lasso_configurations_per_training_size,
                                        job_csv_path)
        else:
            perform_spr_pipeline(training_size_options, brlen_generators, curr_msa_stats,
                                 lasso_configurations_per_training_size,
                                 job_csv_path)
    with open(curr_job_status_file, 'w') as job_status_f:
        job_status_f.write("Done")


if __name__ == "__main__":
    main()
