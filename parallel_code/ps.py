from help_functions import *
from config import *
import time
import shutil
import sklearn


def generate_results_folder(curr_run_prefix):
    create_dir_if_not_exists(RESULTS_FOLDER)
    curr_run_prefix = os.path.join(RESULTS_FOLDER, curr_run_prefix)
    create_or_clean_dir(curr_run_prefix)
    return curr_run_prefix


def distribute_MSAs_over_jobs(file_path_list, all_jobs_results_folder, args):
    jobs_csv_path_list = []
    status_file_path_list = []
    files_per_job = int(len(file_path_list) / args.n_jobs)
    logging.info(
        "Distributing MSAs over jobs: number of files per job={files_per_job} ".format(files_per_job=files_per_job))
    for job_ind in range(min(args.n_jobs, len(file_path_list))):
        curr_job_folder = os.path.join(all_jobs_results_folder, "job_" + str(job_ind))
        create_or_clean_dir(curr_job_folder)
        first_msa_ind = files_per_job * job_ind
        last_msa_ind = files_per_job * (job_ind + 1)
        job_related_files_paths = get_job_related_files_paths(curr_job_folder, job_ind)
        job_msa_paths = file_path_list[first_msa_ind:last_msa_ind]
        jobs_csv_path_list.append(job_related_files_paths["job_csv_path"])
        status_file_path_list.append(job_related_files_paths["job_status_file"])
        raw_results = pd.DataFrame(
        )
        raw_results.to_csv(job_related_files_paths["job_csv_path"], index=False)
        with open(job_related_files_paths["job_msa_paths_file"], 'w') as f:
            for path in job_msa_paths:
                f.write("%s\n" % path)
        logging.info("job number {} will run on files {}".format(job_ind, job_msa_paths))

        run_command = f' python {MAIN_CODE_PATH} ' \
                ' --job_ind {job_ind} --curr_job_folder {curr_job_folder} {previous_args}' \
                .format(
                job_ind=job_ind, previous_args=generate_argument_str(args), curr_job_folder=curr_job_folder
                )
        job_name = args.jobs_prefix + str(job_ind)
        if not LOCAL_RUN:
            submit_linux_job(job_name, curr_job_folder, run_command, 1, 1, job_ind)
        else:
            submit_local_job(MAIN_CODE_PATH, ["--job_ind", str(job_ind), "--curr_job_folder", curr_job_folder
             ] + generate_argument_list(args))
    csv_path_to_status_path_dict = {csv_path: status_path for csv_path, status_path in
                                    zip(jobs_csv_path_list, status_file_path_list)}
    return csv_path_to_status_path_dict


def main():

    parser = main_parser()
    args = parser.parse_args()
    generate_argument_str(args)
    all_jobs_results_folder = generate_results_folder(args.run_prefix)
    all_jobs_general_log_file = os.path.join(all_jobs_results_folder, "log_file.log")
    logging.basicConfig(filename=all_jobs_general_log_file, level=LOGGING_LEVEL)
    logging.info("Args = {args}".format(args=args))
    all_jobs_csv = os.path.join(all_jobs_results_folder, OUTPUT_CSV_NAME + '.csv')
    all_jobs_backup_csv = os.path.join(all_jobs_results_folder, "backup.csv")
    logging.info('#Started running')
    if args.alternative_analysis:
        file_path_list = extract_alignment_files_from_dir(args.alternative_files_folder)
    else:
        file_path_list = extract_alignment_files_from_general_csv(MSAs_CSV_PATH)
    logging.info("There are overall {nMSAs} available ".format(nMSAs=len(file_path_list)))
    if os.path.exists(all_jobs_backup_csv) and os.path.os.stat(all_jobs_backup_csv).st_size > 0:
        shutil.copy(all_jobs_backup_csv, all_jobs_csv)
        file_path_list = [f for f in file_path_list if f not in pd.read_csv(all_jobs_backup_csv)["dataset_id"].unique()]
        logging.info(
            "After removing files that exist in {} there are {} MSAs".format(all_jobs_backup_csv, len(file_path_list)))
    file_path_list = remove_MSAs_with_not_enough_seq(file_path_list, args.min_n_seq)
    logging.info("There are {} MSAs with at least {} sequences".format(len(file_path_list), args.min_n_seq))
    file_path_list = file_path_list[args.first_msa_ind:args.first_msa_ind + args.n_MSAs]
    logging.debug("Alignment files are " + str(file_path_list))
    csv_path_to_status_path_dict = distribute_MSAs_over_jobs(file_path_list, all_jobs_results_folder, args)
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
        time.sleep(WAITING_TIME_UPDATE)
    remove_empty_columns(all_jobs_csv)
    remove_empty_columns(all_jobs_backup_csv)
    logging.info(f'All csvs are updated to {all_jobs_csv}. Done!!!')


if __name__ == "__main__":
    main()
