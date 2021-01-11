
from help_functions import *

#unified_csv_path = '/groups/pupko/noaeker/lasso_positions_sampling_results/current_spr.csv'
#job_folder = '/groups/pupko/noaeker/lasso_positions_sampling_results/brlen_spr'
# log_file_path = '/groups/pupko/noaeker/lasso_positions_sampling_results/add_csvs_log'
from shutil import copyfile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--unified_csv_path', action='store', type=str)
    parser.add_argument('--jobs_folder', action='store', type=str)
    parser.add_argument('--n_jobs', action='store', type=int)
    parser.add_argument('--log_file_path', action='store', type=str)
    args = parser.parse_args()
    logging.basicConfig(filename=args.log_file_path, level=LOGGING_LEVEL)
    pd.DataFrame().to_csv(args.unified_csv_path)
    csv_paths = [os.path.join(args.jobs_folder, "job_{}".format(i), "{}.csv".format(i)) for i in range(args.n_jobs)]
    for csv_path in csv_paths :
        logging.info("csv_path={}".format(csv_path))
        logging.info("path exists: {} and size is : {} can read: {}".format(os.path.exists(csv_path), os.path.getsize(csv_path), os.access(csv_path, os.R_OK)))
    logging.info("copied csv paths are: {}".format(csv_paths ))
    add_csvs_content(csv_paths , args.unified_csv_path)

if __name__ == "__main__":
    main()