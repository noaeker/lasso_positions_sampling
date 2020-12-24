from help_functions import *
from config import *
import subprocess
import sys
import time
import shutil


def generate_results_folder(curr_run_prefix):
    create_dir_if_not_exists(RESULTS_FOLDER)
    curr_run_prefix = os.path.join(RESULTS_FOLDER, curr_run_prefix)
    create_or_clean_dir(curr_run_prefix)
    return curr_run_prefix




def distribute_MSAs_over_jobs(path_list, n_jobs, all_jobs_results_folder, max_n_sequences,
                              n_random_starting_trees, jobs_prefix,only_evaluate_lasso):
    jobs_csv_path_list = []
    status_file_path_list = []
    files_per_job = int(len(path_list) / n_jobs)
    logging.info(
        "Distributing MSAs over jobs: number of files per job={files_per_job} ".format(files_per_job=files_per_job))
    for job_ind in range(min(n_jobs, len(path_list))):
        curr_job_folder = os.path.join(all_jobs_results_folder, "job_" + str(job_ind))
        create_or_clean_dir(curr_job_folder)
        first_msa_ind = files_per_job * job_ind
        last_msa_ind = files_per_job * (job_ind + 1)
        job_msa_paths = path_list[first_msa_ind:last_msa_ind]
        job_msa_paths_file = os.path.join(curr_job_folder, "file_paths_" + str(job_ind))
        job_csv_path = os.path.join(curr_job_folder, str(job_ind) + ".csv")
        jobs_csv_path_list.append(job_csv_path)
        job_status_file = os.path.join(curr_job_folder, str(job_ind) + "_status")
        status_file_path_list.append(job_status_file)
        raw_results = pd.DataFrame(
        )
        raw_results.to_csv(job_csv_path, index=False)
        spr_log_path = os.path.join(curr_job_folder, "job_" + str(job_ind) + "_spr_log.log")
        general_log_path = os.path.join(curr_job_folder, "job_" + str(job_ind) + "_general_log.log")
        with open(job_msa_paths_file, 'w') as f:
            for path in job_msa_paths:
                f.write("%s\n" % path)
        logging.info("job number {} will run on files {}".format(job_ind, job_msa_paths))
        cmds_path = os.path.join(curr_job_folder, str(job_ind) + ".cmds")
        job_log_path = os.path.join(curr_job_folder, str(job_ind) + "_tmp_log")
        job_line = f'module load gcc/gcc-8.2.0; module load python/python-anaconda3.6.5-orenavr2!@#python; python /groups/pupko/noaeker/lasso_positions_sampling/parallel_code/MSA_positions_sampling.py {job_ind} {curr_job_folder} {job_msa_paths_file} {job_csv_path} {spr_log_path} {general_log_path} {job_status_file} {max_n_sequences} {n_random_starting_trees} {only_evaluate_lasso}\t{jobs_prefix}{str(job_ind)}'
        with open(cmds_path, 'w') as cmds_f:
            cmds_f.write(job_line)
        if not LOCAL_RUN:
            os.system(f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {cmds_path} {job_log_path}')
        else:
            msa_code_location = MAIN_CODE_PATH
            theproc = subprocess.Popen(
                [sys.executable, msa_code_location, str(job_ind), curr_job_folder,
                 job_msa_paths_file, job_csv_path, spr_log_path, general_log_path,
                 job_status_file, str(max_n_sequences), str(n_random_starting_trees),str(only_evaluate_lasso)])
            theproc.communicate()
    csv_path_to_status_path_dict = {csv_path: status_path for csv_path, status_path in
                                    zip(jobs_csv_path_list, status_file_path_list)}
    return csv_path_to_status_path_dict



def main():
    args = ps_parser()
    all_jobs_results_folder = generate_results_folder(args.run_prefix)
    all_jobs_general_log_file = os.path.join(all_jobs_results_folder, "log_file.log")
    logging.basicConfig(filename=all_jobs_general_log_file, level=LOGGING_LEVEL)
    logging.info(
        "Program input:\nn_MSAs = {}\nn_jobs = {}\nMSA_start_ind = {}\nmax_n_sequences = {}\nn_random_starting_tree = {}\nonly_evaluate_lasso={}".format(
            args.n_MSAs, args.n_jobs, args.first_msa_ind, args.max_n_seq, args.n_random_starting_trees,
            args.only_evaluate_lasso))
    all_jobs_csv = os.path.join(all_jobs_results_folder, OUTPUT_CSV_NAME + '.csv')
    all_jobs_backup_csv = os.path.join(all_jobs_results_folder, "backup.csv")
    logging.info('#Started running')
    file_path_list = extract_alignment_files_from_general_csv(MSAs_CSV_PATH)
    logging.info("There are overall {nMSAs} available ".format(nMSAs=len(file_path_list)))
    if os.path.exists(all_jobs_backup_csv) and os.path.os.stat(all_jobs_backup_csv).st_size > 0:
        shutil.copy(all_jobs_backup_csv, all_jobs_csv)
        file_path_list = [f for f in file_path_list if f not in pd.read_csv(all_jobs_backup_csv)["dataset_id"].unique()]
        logging.info(
            "After removing files that exist in {} there are {} MSAs".format(all_jobs_backup_csv, len(file_path_list)))
    file_path_list = remove_MSAs_with_not_enough_seq(file_path_list, MIN_N_SEQ)
    logging.info("There are {} MSAs with at least {} sequences".format(len(file_path_list), MIN_N_SEQ))
    file_path_list = file_path_list[args.first_msa_ind:args.first_msa_ind + args.n_MSAs]
    logging.debug("Alignment files are " + str(file_path_list))
    csv_path_to_status_path_dict = distribute_MSAs_over_jobs(file_path_list, args.n_jobs,
                                                             all_jobs_results_folder, args.max_n_seq,
                                                             args.n_random_starting_trees, args.jobs_prefix,args.only_evaluate_lasso)
    while len(csv_path_to_status_path_dict) > 0:
        # logging.info(f'***Checking CSVs status***')
        new_csvs_to_update = []
        for csv_path in csv_path_to_status_path_dict:  # check which csvs are ready to be updated
            status_path = csv_path_to_status_path_dict[csv_path]
            if os.path.exists(status_path):
                new_csvs_to_update.append(csv_path)
        for key in new_csvs_to_update:
            if key in csv_path_to_status_path_dict:
                del csv_path_to_status_path_dict[key]
        if len(new_csvs_to_update) > 0:
            t = time.localtime()
            current_time = time.strftime("%m/%d/%Y, %H:%M:%S", t)
            logging.info(
                "Current time: {} About to update following jobs local csvs to general csv: {}".format(current_time,
                                                                                                       new_csvs_to_update))
            logging.info("Number of jobs still running is : {} ".format(len(csv_path_to_status_path_dict)))
            add_csvs_content(new_csvs_to_update, all_jobs_csv)
            add_csvs_content(new_csvs_to_update, all_jobs_backup_csv)
        time.sleep(WAITING_TIME_CSV_UPDATE)
    remove_empty_columns(all_jobs_csv)
    remove_empty_columns(all_jobs_backup_csv)
    logging.info(f'All csvs are updated to {all_jobs_csv}. Done!!!')


if __name__ == "__main__":
    main()
