
from generate_SPR import *
from training_and_test_set_generation import *
from raxml import *
import pickle


def generate_or_copy_random_starting_tree(i, curr_run_directory, curr_msa_stats):
    seed=SEED+i
    random_tree_folder = os.path.join(curr_run_directory, "RANDOM_starting_tree_" + str(i))
    create_or_clean_dir(random_tree_folder)
    starting_tree_path = os.path.join(random_tree_folder, ".raxml.startTree")
    baseline_starting_tree_path = starting_tree_path.replace(curr_msa_stats["run_prefix"],
                                                             curr_msa_stats["spr_baseline_run_prefix"])
    if os.path.exists(baseline_starting_tree_path):
        logging.info("Using starting tree in {}".format(baseline_starting_tree_path))
        shutil.copyfile(baseline_starting_tree_path, starting_tree_path)
    if not os.path.exists(starting_tree_path):
        logging.info("Generating a totally random tree as a starting tree")
        starting_tree_path,elapsed_running_time = generate_n_random_tree_topology_constant_brlen(curr_msa_stats = curr_msa_stats,n=1, alpha=curr_msa_stats["alpha"],
                                                                            original_file_path=curr_msa_stats["local_alignment_path"],
                                                                            curr_run_directory=random_tree_folder, seed=seed)
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
    constant_sites_pct = len((alignment_df.loc[:,alignment_df.apply(pd.Series.nunique) == 1]).columns)/len(alignment_df.columns)
    probabilities = [list(map(lambda x: x / n_seq, counts_per_position[col].values())) for col in
                     list(alignment_df)]
    entropy = [sum(list(map(lambda x: -x * np.log(x), probabilities[col]))) for col in list(alignment_df)]
    avg_entropy = np.mean(entropy)
    return constant_sites_pct, avg_entropy, gap_positions_pct


def generate_msa_general_stats( original_alignment_path, file_ind, curr_msa_version_folder,args
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
    local_full_msa_path = os.path.join(curr_msa_version_folder, file_name + file_type)
    take_up_to_x_sequences(original_alignment_data, local_full_msa_path, args.max_n_seq, file_type_biopython, args.max_n_loci)
    with open(local_full_msa_path) as original:
        reduced_local_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    reduced_local_alignment_df = alignment_list_to_df(reduced_local_alignment_data)
    n_seq, n_loci = reduced_local_alignment_df.shape
    constant_sites_pct,avg_entropy, gap_positions_pct = get_positions_stats(reduced_local_alignment_df, n_seq)
    curr_msa_stats = {"n_seq": n_seq, "MSA_original_n_seq": orig_n_seq,
                      "n_seq_before_reduction_by_RaxML": n_seq,
                      "n_loci": n_loci, "MSA_original_n_loci": n_loci, "dataset_id": dataset_id,
                      "curr_msa_version_folder": curr_msa_version_folder,
                      "alignment_data": reduced_local_alignment_data,
                      "MSA_original_alignment_data": original_alignment_data,
                      "local_alignment_path": local_full_msa_path,
                      "full_data_unique_name": full_data_unique_name,
                      "file_ind": file_ind,
                      "file_name": file_name,
                      "file_type": file_type,
                      "file_type_biopython": file_type_biopython,
                      "constant_sites_pct": constant_sites_pct,
                      "avg_entropy": avg_entropy,
                      "gap_pct": gap_positions_pct,


                      }
    curr_msa_stats.update(vars(args))
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
    logging.info("Reduced MSA stats computed:\n {reduced_stats}".format(
        reduced_stats={key: reduced_curr_msa_stats[key] for key in reduced_curr_msa_stats if key not in IGNORE_COLS_IN_CSV}))

    return reduced_curr_msa_stats


def update_chosen_brlen_generators(exp_brlen, uni_brlen, opt_brlen, const_brlen):
    brlen_generators = {}
    for brlen_generator_name, brlen_generator_func in BRLEN_GENERATORS.items():
        if (brlen_generator_name == 'exponential' and exp_brlen) or (
                brlen_generator_name == 'uniform' and uni_brlen) or (
                brlen_generator_name == 'optimized' and opt_brlen) or (
                brlen_generator_name == 'const' and const_brlen):
            brlen_generators[brlen_generator_name] = brlen_generator_func
    return brlen_generators
    logging.info("Brlen genertors are chosen to be {}".format(str(brlen_generators.keys())))






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
            with open(lasso_results["sampled_alignment_path"]) as sampled_path:
                sampled_data = list(SeqIO.parse(sampled_path, curr_msa_stats["file_type_biopython"]))
            sampled_alignment_df = alignment_list_to_df( sampled_data )
            n_seq, n_loci = sampled_alignment_df.shape
            constant_sites_pct, avg_entropy, gap_positions_pct = get_positions_stats(sampled_alignment_df, n_seq)
            results_dict = {"constant_sites_pct_sampled" : constant_sites_pct,"avg_entropy_sampled" : avg_entropy, "gap_positions_pct_sampled":gap_positions_pct }
            logging.info(f"Sampled data metrics: {results_dict}")
            curr_msa_stats.update(results_dict)
            curr_msa_stats.update(lasso_results)
            logging.info("only evaluating lasso on brlen {} and training size {}: ".format(brlen_generator_name, training_size))
            lasso_evaluation_result = {k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                       k not in IGNORE_COLS_IN_CSV
                                       }
            all_msa_results = all_msa_results.append(lasso_evaluation_result, ignore_index=True)
            all_msa_results.to_csv(job_csv_path, index=False)


def dilute_msa(curr_msa_stats, curr_run_directory):
    dilute_amount = curr_msa_stats["dilute_amount"]
    sampled_alignment_path_extended = os.path.join(curr_run_directory,"diluted_sampled_msa")
    weights_path_extended = os.path.join(curr_run_directory,"diluted_weights")
    logging.info('Diluting MSA with {} positions, writing results to {} and writing weights to {}'.format(dilute_amount,
                                                                                                          sampled_alignment_path_extended,
                                                                                                          weights_path_extended))
    total_lasso_weights = sum(curr_msa_stats["lasso_chosen_weights"])
    random.seed(SEED)
    diluted_samp_indexes = curr_msa_stats["lasso_chosen_locis"] + random.sample(list(range(curr_msa_stats["n_loci"])),
                                                                                curr_msa_stats["dilute_amount"])
    write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], sampled_alignment_path_extended,
                                    diluted_samp_indexes, curr_msa_stats["file_type_biopython"])
    curr_msa_stats["sampled_alignment_path"] = sampled_alignment_path_extended
    extra_weights = [total_lasso_weights / (dilute_amount * curr_msa_stats["dilute_mul"])] * dilute_amount
    with open(weights_path_extended, 'w') as f:
        for weight in curr_msa_stats["lasso_chosen_weights"] + extra_weights:
            if USE_INTEGER_WEIGHTS:
                weight= int(weight)
            f.write(str(weight) + " ")
    curr_msa_stats["weights_file_path"] = weights_path_extended

def perform_raxml_search_pipeline(training_size_options, brlen_generators, curr_msa_stats, lasso_configurations_per_training_size,
                                job_csv_path):
    all_msa_results = pd.DataFrame(
    )
    all_msa_results.to_csv(job_csv_path, index=False)
    curr_run_directory = os.path.join(curr_msa_stats["curr_msa_version_folder"],"RaxML_search")
    create_dir_if_not_exists(curr_run_directory)
    standard_run_folder= os.path.join(curr_run_directory,"standard_run")
    create_dir_if_not_exists(standard_run_folder)
    standard_raxml_results_dump = os.path.join(standard_run_folder, 'standard_RAxML.dump')
    standard_raxml_dump_baseline = standard_raxml_results_dump.replace(curr_msa_stats["run_prefix"],
                                                                     curr_msa_stats["RAxML_baseline_run_prefix"])
    if os.path.exists(standard_raxml_dump_baseline):
        logging.info("Using dump standard RAxML results in {}".format(standard_raxml_dump_baseline))
        with open(standard_raxml_dump_baseline, 'rb') as handle:
            standard_raxml_search_results = pickle.load(handle)
    else:
        logging.info("Performing standard RAxML search and saving results to {}".format(standard_raxml_results_dump))
        standard_raxml_search_results = raxml_search_pipeline(standard_run_folder, curr_msa_stats,
                                                              curr_msa_stats["n_raxml_parsimony_trees"], curr_msa_stats["n_raxml_random_trees"],
                                                              standrad_search=True)
        with open(standard_raxml_results_dump, 'wb') as handle:
            pickle.dump(standard_raxml_search_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
    curr_msa_stats.update(standard_raxml_search_results)
    logging.info("Standard RAxML results: \n{}".format(standard_raxml_search_results))
    if not curr_msa_stats["only_full_search"]:
        logging.info("**Now performing a Lasso-based search")
        for brlen_generator_name in brlen_generators:
            brlen_run_directory = os.path.join(curr_run_directory, brlen_generator_name)
            create_dir_if_not_exists(brlen_run_directory)
            curr_msa_stats["brlen_generator"] = brlen_generator_name
            for training_size in training_size_options:
                curr_msa_stats["actucal_training_size"] = training_size
                curr_training_size_and_brlen_directory = os.path.join(brlen_run_directory, str(training_size))
                create_dir_if_not_exists(curr_training_size_and_brlen_directory)
                lasso_thresholds_during_search = [float(t) for t in curr_msa_stats['lasso_thresholds_search'].split("_")]
                for threshold in lasso_thresholds_during_search:
                    curr_threshold_directory = os.path.join(curr_training_size_and_brlen_directory , str(threshold))
                    create_dir_if_not_exists(curr_threshold_directory)
                    lasso_config_per_threshold = lasso_configurations_per_training_size[brlen_generator_name][
                        training_size]
                    if threshold in lasso_config_per_threshold:
                        lasso_results = lasso_config_per_threshold[threshold]
                    else:
                        continue
                    curr_msa_stats.update(lasso_results)
                    logging.info("Starting Lasso-based RaxML search using {brlen} brlen and training size: {size}".format(size=training_size, brlen = brlen_generator_name ))
                    lasso_based_RAxML_results = raxml_search_pipeline(curr_threshold_directory,curr_msa_stats, curr_msa_stats["n_raxml_parsimony_trees"], curr_msa_stats["n_raxml_random_trees"], standrad_search= False)
                    logging.info(f'lasso based RAxML results are:\n{lasso_based_RAxML_results}')
                    curr_msa_stats.update(lasso_based_RAxML_results)
                    best_tree_full_data_path = curr_msa_stats["standard_best_tree_path"]
                    best_tree_first_phase = curr_msa_stats["lasso_first_phase_best_tree"]
                    rf_folder = os.path.join(curr_threshold_directory,"rf_calculations")
                    create_dir_if_not_exists(rf_folder)
                    rf_first_phase_trees =unify_text_files([best_tree_full_data_path,best_tree_first_phase],os.path.join(rf_folder,"rf_first_phase_trees"))
                    curr_msa_stats["rf_dist_first_phase"] = calculate_rf_dist(rf_first_phase_trees, rf_folder,prefix="rf_first_phase")
                    if curr_msa_stats["do_raxml_lasso_second_phase"]:
                        best_tree_second_phase = curr_msa_stats["lasso_second_phase_best_tree"]
                        rf_second_phase_trees = unify_text_files([best_tree_full_data_path, best_tree_second_phase],os.path.join(rf_folder,"rf_second_phase_trees"))
                        curr_msa_stats["rf_dist_second_phase"] = calculate_rf_dist(rf_second_phase_trees, rf_folder,prefix="rf_second_phase")
                    all_msa_results = all_msa_results.append({k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                                              k not in IGNORE_COLS_IN_CSV
                                                              }, ignore_index=True)

    else:
        logging.info("Updating full MSA results to csv and finishing")
        all_msa_results = all_msa_results.append({k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                                  k not in IGNORE_COLS_IN_CSV
                                                  }, ignore_index=True)
    all_msa_results.to_csv(job_csv_path)

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
            "Generating folder for starting tree {i} in dicrectory {dir}".format(i=i, dir=starting_tree_run_directory))
        starting_tree_path = generate_or_copy_random_starting_tree(i, starting_tree_run_directory, curr_msa_stats)
        curr_msa_stats["starting_tree_path"] = starting_tree_path
        naive_spr_results_dump = os.path.join(starting_tree_run_directory, 'naive_spr.dump')
        naive_spr_results_dump_baseline = naive_spr_results_dump.replace(curr_msa_stats["run_prefix"],curr_msa_stats["spr_baseline_run_prefix"])
        if os.path.exists(naive_spr_results_dump_baseline):
            logging.info("Using naive SPR results in {}".format(naive_spr_results_dump_baseline))
            with open(naive_spr_results_dump_baseline, 'rb') as handle:
                naive_spr_results = pickle.load(handle)
        else:
            naive_spr_results = naive_spr_on_current_starting_tree(starting_tree_path, curr_msa_stats,
                                                                   curr_run_directory=starting_tree_run_directory)
            with open(naive_spr_results_dump, 'wb') as handle:
                pickle.dump(naive_spr_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        curr_msa_stats.update(naive_spr_results)
        lasso_thresholds_during_search = [float(t) for t in curr_msa_stats['lasso_thresholds_search'].split("_")]
        for brlen_generator_name in brlen_generators:
            brlen_run_directory = os.path.join(starting_tree_run_directory, brlen_generator_name)
            create_dir_if_not_exists(brlen_run_directory)
            curr_msa_stats["brlen_generator"] = brlen_generator_name
            for training_size in training_size_options:
                curr_msa_stats["actucal_training_size"] = training_size
                curr_training_size_directory = os.path.join(brlen_run_directory, str(training_size))
                create_dir_if_not_exists(curr_training_size_directory)
                for threshold in  lasso_thresholds_during_search:
                    curr_threshold_directory = os.path.join(curr_training_size_directory , str(threshold))
                    create_dir_if_not_exists( curr_threshold_directory )
                    lasso_config_per_threshold = lasso_configurations_per_training_size[brlen_generator_name][training_size]
                    if threshold in lasso_config_per_threshold:
                        lasso_results = lasso_config_per_threshold[threshold]
                    else:
                        continue
                    curr_msa_stats.update(lasso_results)
                    logging.info("Starting Lasso-based SPR on training size: {size}".format(size=training_size))
                    start_time = time.time()
                    lasso_based_spr_results = lasso_based_spr_on_current_starting_tree(starting_tree_path,
                                                                                       curr_msa_stats,
                                                                                       curr_threshold_directory)
                    spr_pipeline_running_time = time.time()-start_time
                    curr_msa_stats.update(lasso_based_spr_results)

                    best_tree_full_newick= curr_msa_stats["naive_SPR_tree_newick"]
                    best_tree_first_newick= curr_msa_stats["lasso_SPR_first_phase_tree_newick" ]
                    best_tree_second_phase_newick = curr_msa_stats["lasso_SPR_second_phase_tree_newick"]
                    rf_folder = os.path.join( curr_threshold_directory, "rf_calculations")
                    create_dir_if_not_exists(rf_folder)
                    rf_first_phase_trees = os.path.join(rf_folder, "rf_first_phase_trees")
                    with open(rf_first_phase_trees,'w') as FIRST_PHASE_RF:
                        FIRST_PHASE_RF.writelines([best_tree_full_newick,"\n", best_tree_first_newick])
                    rf_second_phase_trees = os.path.join(rf_folder, "rf_second_phase_trees")
                    with open( rf_second_phase_trees, 'w') as SECOND_PHASE_RF:
                        SECOND_PHASE_RF.writelines([best_tree_full_newick,"\n", best_tree_second_phase_newick ])

                    curr_msa_stats["rf_dist_first_phase"] = calculate_rf_dist(rf_first_phase_trees, rf_folder, prefix = "rf_first_phase")
                    curr_msa_stats["rf_dist_second_phase"] = calculate_rf_dist(rf_second_phase_trees, rf_folder, prefix="rf_second_phase")

                    curr_msa_stats.update({'spr_pipeline_running_time':spr_pipeline_running_time})

                    all_msa_results = all_msa_results.append({k: curr_msa_stats[k] for k in curr_msa_stats.keys() if
                                                              k not in IGNORE_COLS_IN_CSV
                                                              }, ignore_index=True)
                    all_msa_results.to_csv(job_csv_path)







def get_lasso_configurations(curr_msa_version_folder,args,brlen_generators, curr_msa_stats,training_size_options):
    curr_msa_version_lasso_dump = os.path.join(curr_msa_version_folder, 'lasso.dump')
    curr_msa_version_lasso_dump_baseline = curr_msa_version_lasso_dump.replace(args.run_prefix,
                                                                               args.lasso_baseline_run_prefix)
    if os.path.exists(curr_msa_version_lasso_dump_baseline):
        with open(curr_msa_version_lasso_dump_baseline, 'rb') as handle:
            logging.info(
                "Using lasso dump files in {} ".format(curr_msa_version_lasso_dump))
            lasso_configurations_per_training_size = pickle.load(handle)
    else:
        lasso_configurations_per_training_size = Lasso_training_and_test(brlen_generators, curr_msa_stats,
                                                                         training_size_options,
                                                                         args.random_trees_test_size)
        with open(curr_msa_version_lasso_dump, 'wb') as handle:
            pickle.dump(lasso_configurations_per_training_size, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return lasso_configurations_per_training_size

def get_msa_stats(curr_msa_version_folder, original_alignment_path, args, file_ind):
    curr_msa_version_stats_dump = os.path.join(curr_msa_version_folder, 'curr_msa_stats.dump')
    curr_msa_version_stats_dump_baseline = curr_msa_version_stats_dump.replace(args.run_prefix,
                                                                               args.msa_baseline_run_prefix)

    if os.path.exists(curr_msa_version_stats_dump_baseline):
        with open(curr_msa_version_stats_dump_baseline, 'rb') as handle:
            logging.info(
                "Using msa stats dump files in {} ".format(curr_msa_version_stats_dump_baseline ))
            curr_msa_stats = pickle.load(handle)
            curr_msa_stats["curr_msa_version_folder"] = curr_msa_version_folder
            curr_msa_stats.update(vars(args))
    else:
        logging.info(
            "Generating msa stats from beggining".format(curr_msa_version_stats_dump))
        curr_msa_stats = generate_msa_general_stats(
            original_alignment_path, file_ind, curr_msa_version_folder, args)
        try:
            extract_and_update_RaxML_statistics_from_full_data(
                curr_msa_stats)
        except RE_RUN_ON_REDUCED_VERSION:
            curr_msa_stats = re_run_on_reduced_version(curr_msa_stats, original_alignment_path,
                                                       file_ind
                                                       )  # use current curr_msa_stats
            extract_and_update_RaxML_statistics_from_full_data(
                curr_msa_stats)
        with open(curr_msa_version_stats_dump, 'wb') as handle:
            pickle.dump(curr_msa_stats, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return curr_msa_stats




def main():
    parser = job_parser()
    args = parser.parse_args()
    job_related_file_paths = get_job_related_files_paths(args.curr_job_folder, args.job_ind)
    job_msa_paths_file, general_log_path, job_csv_path, curr_job_status_file = job_related_file_paths[
                                                                                                 "job_msa_paths_file"], \
                                                                                             job_related_file_paths[
                                                                                                 "general_log_path"], \
                                                                                             job_related_file_paths[
                                                                                                 "job_csv_path"], \
                                                                                             job_related_file_paths[
                                                                                                 "job_status_file"]
    with open(job_msa_paths_file, "r") as paths_file:
        curr_job_file_path_list = paths_file.read().splitlines()
    logging.basicConfig(filename=general_log_path, level=LOGGING_LEVEL)
    logging.info('#Started running on job' + str(args.job_ind))
    logging.info("Job arguments : {}".format(args))
    training_size_options = [int(size) for size in args.random_trees_training_size.split("_")]
    brlen_generators = update_chosen_brlen_generators(args.exp_brlen, args.uni_brlen, args.opt_brlen, args.const_brlen)
    for file_ind, original_alignment_path in enumerate(curr_job_file_path_list):
        msa_name = original_alignment_path.replace(MSAs_FOLDER, "").replace("ref_msa.aa.phy", "").replace(os.path.sep,
                                                                                                          "_")
        logging.info('#running on file name {} and ind (relativ to job) {} path= {}'.format(msa_name, file_ind,
                                                                                            original_alignment_path))
        curr_msa_version_folder = os.path.join(args.curr_job_folder, msa_name)
        create_or_clean_dir(curr_msa_version_folder)
        curr_msa_stats = get_msa_stats(curr_msa_version_folder, original_alignment_path, args, file_ind)
        if not args.only_full_search:
            lasso_configurations_per_training_size = get_lasso_configurations(curr_msa_version_folder,args,brlen_generators, curr_msa_stats,training_size_options)
        else:
            logging.info("Not calculating or using existing Lasso results")
            lasso_configurations_per_training_size = None
        if args.only_evaluate_lasso:
            perform_only_lasso_pipeline(training_size_options, brlen_generators, curr_msa_stats,
                                        lasso_configurations_per_training_size,
                                        job_csv_path)
        elif args.use_raxml_search:
            perform_raxml_search_pipeline(training_size_options, brlen_generators, curr_msa_stats,
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
