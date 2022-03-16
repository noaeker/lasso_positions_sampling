
from lasso_evaluation_pipeline import *
from SPR_pipeline import perform_spr_pipeline
from MSA_features_extraction_pipeline import *
from config import BRLEN_GENERATORS, EXAMPLE_RUN
from help_functions import get_job_related_files_paths, create_or_clean_dir
from input_parsers import job_parser
import pandas as pd
from config import LOGGING_LEVEL,CSV_MSAs_FOLDER



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
    n_seq_options = [int(t) for t in
                     (args.max_n_seq).split("_")]
    n_loci_options = [int(t) for t in
                      (args.max_n_loci).split("_")]
    all_msa_results = pd.DataFrame(
    )
    for file_ind, original_alignment_path in enumerate(curr_job_file_path_list):
        msa_name = os.path.basename(original_alignment_path) if EXAMPLE_RUN else original_alignment_path.replace(CSV_MSAs_FOLDER, "").replace("example_msa.phy", "").replace(os.path.sep,
                                                                                                          "_")
        logging.info(
            f'#running on file name {msa_name} and ind (relativ to job) {file_ind}  original path= {original_alignment_path}')
        curr_msa_folder = os.path.join(args.curr_job_folder, msa_name)
        create_or_clean_dir(curr_msa_folder)
        all_msa_results.to_csv(job_csv_path, index=False, sep ='\t')
        for actual_n_seq in n_seq_options:
            curr_n_seq_folder = os.path.join(curr_msa_folder, f"n_seq_{actual_n_seq}")
            create_or_clean_dir(curr_n_seq_folder)
            for actual_n_loci in n_loci_options:
                    curr_n_loci_folder = os.path.join(curr_n_seq_folder, f"n_loci_{actual_n_loci}")
                    create_or_clean_dir(curr_n_loci_folder)
                    logging.info(f" Trimming MSA to n_seq = {actual_n_seq} n_loci = {actual_n_loci}")
                    curr_msa_stats = get_msa_stats(curr_n_loci_folder, original_alignment_path, args, file_ind,
                                                   actual_n_seq, actual_n_loci)
                    if curr_msa_stats == -1:
                        continue
                    curr_msa_stats["actual_n_loci"] = actual_n_loci
                    curr_msa_stats["actual_n_seq"] = actual_n_seq
                    if not args.only_full_search:
                        lasso_configurations_per_training_size = get_lasso_configurations(curr_n_loci_folder, args,
                                                                                          brlen_generators,
                                                                                          curr_msa_stats,
                                                                                          training_size_options)
                    else:
                        lasso_configurations_per_training_size = None
                    if args.only_evaluate_lasso:
                        all_msa_results = perform_only_lasso_pipeline(training_size_options, brlen_generators,
                                                                      curr_msa_stats,
                                                                      lasso_configurations_per_training_size,
                                                                      job_csv_path, all_msa_results, curr_n_loci_folder)
                    else:
                        all_msa_results = perform_spr_pipeline(training_size_options, brlen_generators, curr_msa_stats,
                                                               lasso_configurations_per_training_size,
                                                               job_csv_path, all_msa_results)

    with open(curr_job_status_file, 'w') as job_status_f:
        job_status_f.write("Done")


if __name__ == "__main__":
    main()
