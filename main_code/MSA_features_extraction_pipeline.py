
from raxml import extract_raxml_statistics_from_msa,RE_RUN_ON_REDUCED_VERSION
import pickle
import logging
import os
from help_functions import delete_dir_content,extract_file_type, trim_MSA,alignment_list_to_df,get_positions_stats
from Bio import SeqIO
from partitioned_analysis import parse_raxml_partition_file,generate_loci_corrected_partition_model_file
from config import PARTITION_MODELS_FILE,IGNORE_COLS_IN_CSV
import numpy as np
from rate4site import parse_rate4site,get_rate4site


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


def generate_msa_general_stats(original_alignment_path, file_ind, curr_msa_version_folder, args, actual_n_seq,
                               actual_n_loci
                               ):
    dataset_id = original_alignment_path
    file_name = str(file_ind)
    full_data_unique_name = file_name
    file_type = extract_file_type(original_alignment_path, False)
    file_type_biopython = extract_file_type(original_alignment_path, True)
    with open(original_alignment_path) as original:
        original_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    orig_n_seq = len(original_alignment_data)
    orig_n_loci = len(original_alignment_data[0])
    local_full_msa_path = os.path.join(curr_msa_version_folder, file_name + file_type)
    per_loci_partition = None
    msa_short_name = original_alignment_path.split('_')[-1].split(".")[0]
    msa_model_file = os.path.join(PARTITION_MODELS_FILE, f"{msa_short_name}.raxml.model")
    if os.path.exists(msa_model_file) and (args.compare_loci_gene_distribution or args.do_partitioned_lasso_analysis):
            per_loci_partition,partition_ind_to_name = parse_raxml_partition_file(msa_model_file,  orig_n_loci)
    if orig_n_seq >= actual_n_seq or orig_n_loci>=actual_n_loci:
        logging.info(f"Original number of sequences is {orig_n_seq} and it will be trimmed to {actual_n_seq}\nOriginal number of loci's' is {orig_n_loci} and it will be trimmed to {actual_n_loci}")
        corrected_partition_results = trim_MSA(original_alignment_data, local_full_msa_path, actual_n_seq, file_type_biopython,
                 actual_n_loci, args.loci_shift, per_loci_partition)
        corrected_partition_file = generate_loci_corrected_partition_model_file(corrected_partition_results,partition_ind_to_name, curr_run_directory = curr_msa_version_folder)
    else:
        return -1
    with open(local_full_msa_path) as original:
        reduced_local_alignment_data = list(SeqIO.parse(original, file_type_biopython))
    reduced_local_alignment_df = alignment_list_to_df(reduced_local_alignment_data)
    n_seq, n_loci = reduced_local_alignment_df.shape
    constant_sites_pct, avg_entropy, gap_positions_pct = get_positions_stats(reduced_local_alignment_df)
    curr_msa_stats = {"n_seq": n_seq, "MSA_original_n_seq": orig_n_seq,
                      "n_seq_before_reduction_by_RaxML": n_seq,
                      "n_loci": n_loci, "MSA_original_n_loci": orig_n_loci, "dataset_id": dataset_id,
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
                      "per_loci_partition": corrected_partition_results,
                      "msa_orig_model_file" : msa_model_file,
                      "msa_corrected_model_file": corrected_partition_file

                      }
    curr_msa_stats.update(vars(args))
    logging.info("Basic MSA stats computed:\n {curr_msa_stats}".format(
        curr_msa_stats={key: curr_msa_stats[key] for key in curr_msa_stats if key not in IGNORE_COLS_IN_CSV}))
    return curr_msa_stats


def re_run_on_reduced_version(curr_msa_stats, file_ind):
    raxml_reduced_file = curr_msa_stats["orig_reduced_file_path"]
    reduced_dir, rediced_fname = os.path.split(raxml_reduced_file)
    raxml_reduced_file_renamed = os.path.join(reduced_dir, curr_msa_stats["file_name"] + "_fixed" + extract_file_type(
        raxml_reduced_file))
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
    constant_sites_pct_reduced, avg_entropy_reduced, gap_pct_reduced = get_positions_stats(
        original_alignment_df_reduced
    )
    reduced_curr_msa_stats = curr_msa_stats.copy()
    update_dict = {"file_name": file_name_reduced, "local_alignment_path": raxml_reduced_file_renamed,
                   "n_seq": n_seq_reduced,
                   "n_loci": n_loci_reduced,
                   "alignment_data": reduced_data,
                   "constant_sites_pct": constant_sites_pct_reduced,
                   "avg_entropy": avg_entropy_reduced,
                   "gap_pct": gap_pct_reduced,
                   "file_type_biopython": file_type_biopython
                   }
    reduced_curr_msa_stats.update(
        update_dict)
    logging.info("Reduced MSA stats computed:\n {reduced_stats}".format(
        reduced_stats={key: reduced_curr_msa_stats[key] for key in reduced_curr_msa_stats if
                       key not in IGNORE_COLS_IN_CSV}))

    return reduced_curr_msa_stats



def get_msa_stats(curr_msa_version_folder, original_alignment_path, args, file_ind, actual_n_seq, actual_n_loci):
    curr_msa_version_stats_dump = os.path.join(curr_msa_version_folder, 'curr_msa_stats.dump')
    curr_msa_version_stats_dump_baseline = curr_msa_version_stats_dump.replace(args.run_prefix,
                                                                               args.msa_baseline_run_prefix)

    if os.path.exists(curr_msa_version_stats_dump_baseline):
        with open(curr_msa_version_stats_dump_baseline, 'rb') as handle:
            logging.info(
                "Using msa stats dump files in {} ".format(curr_msa_version_stats_dump_baseline))
            curr_msa_stats = pickle.load(handle)
            curr_msa_stats["curr_msa_version_folder"] = curr_msa_version_folder
            curr_msa_stats.update(vars(args))
    else:
        logging.info(
            "Generating msa stats from beggining".format(curr_msa_version_stats_dump))
        curr_msa_stats = generate_msa_general_stats(
            original_alignment_path, file_ind, curr_msa_version_folder, args, actual_n_seq, actual_n_loci)
        if curr_msa_stats == -1:
            return -1
        try:
            extract_and_update_RaxML_statistics_from_full_data(
                curr_msa_stats)
        except RE_RUN_ON_REDUCED_VERSION:
            curr_msa_stats = re_run_on_reduced_version(curr_msa_stats,
                                                       file_ind
                                                       )  # use current curr_msa_stats
            extract_and_update_RaxML_statistics_from_full_data(
                curr_msa_stats)
        msa_corrected_model_partition_optimized, partition_ind_to_name_optimized = parse_raxml_partition_file(curr_msa_stats["pars_optimized_model"], actual_n_loci)
        curr_msa_stats.update({"msa_corrected_model_partition_optimized" : msa_corrected_model_partition_optimized, "partition_ind_to_name_optimized": partition_ind_to_name_optimized})
        rate4site_output_path = os.path.join(curr_msa_version_folder, "r4s.res")
        rate4site_command = get_rate4site(curr_msa_stats["local_alignment_path"],
                                          curr_msa_stats["parsimony_optimized_tree_path"], rate4site_output_path)
        if os.path.exists(rate4site_output_path):
            try:
                rate4site_scores = parse_rate4site(rate4site_output_path)
                curr_msa_stats["rate4site_scores"] = rate4site_scores
                curr_msa_stats["mean_rate4site_scores"] = np.mean(rate4site_scores)
                n_rate4site_scores = len(rate4site_scores)
                logging.info(f'Succesfully obtained {n_rate4site_scores} rate4site weights')
            except Exception as e:
                logging.error('Failed to parse rate4site scores: ' + str(e))

        else:
            logging.error("Could not generate rate4site output: command = {command}".format(command=rate4site_command))
        with open(curr_msa_version_stats_dump, 'wb') as handle:
            pickle.dump(curr_msa_stats, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return curr_msa_stats

