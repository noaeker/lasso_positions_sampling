import argparse
from sklearn import linear_model
from user_code.raxml_user_version import *




def apply_lasso_on_data_and_update_stats(curr_msa_stats):
    logging.info('Applying Lasso on sitelh data' )
    logging.info("Computing locis and weight using lasso")
    sitelh_training_df = curr_msa_stats["sitelh_df"]
    y_training = sitelh_training_df.sum(axis=1)
    logging.info("Sitelh df dimensions, for lasso computations, are: " + str(sitelh_training_df.shape))
    lasso_model = linear_model.LassoCV(cv=5, normalize=True, max_iter=100000).fit(sitelh_training_df,
                                                                                  y_training)  # add positive=True if using RaxML
    chosen_locis = [ind for ind in range(len(lasso_model.coef_)) if abs(lasso_model.coef_[ind]) > 0]
    weights_file_path = os.path.join(curr_msa_stats.get("prefix") + (".lasso_weights"))
    logging.info("About to write weights file to : " + weights_file_path)
    with open(weights_file_path, 'w') as f:
        for coef in lasso_model.coef_:
            f.write(str(coef) + " ")
    samp_indexes_pct=len(chosen_locis) / len(sitelh_training_df.columns)
    lasso_results_dict={"number of positions chosen:": len(chosen_locis),"chosen positions percent": samp_indexes_pct}
    logging.info("Lasso statistics: {x} out of {y} positions were chosen by lasso".format(x=len(chosen_locis), y=len(sitelh_training_df.columns)))
    curr_msa_stats.update( lasso_results_dict)





def generate_site_lh_data(curr_msa_stats, n_iter):
    msa_path = curr_msa_stats["msa_path"]
    sitelh_csv_path = os.path.join(curr_msa_stats["prefix"], "sitelh.csv")
    logging.info(
        "Computing sitelh data")
    logging.info("Generating " + str(n_iter) + " random trees ")
    sitelh_ll_list = []
    random_trees_folder = curr_msa_stats["prefix"]
    for i in range(n_iter):
        curr_site_lh = raxml_compute_per_site_ll_on_a_random_tree(msa_path, i, curr_msa_stats,random_trees_folder=random_trees_folder)
        sitelh_ll_list.append(curr_site_lh)
    sitelh_df = pd.DataFrame(sitelh_ll_list, columns=list(range(len(sitelh_ll_list[0]))),
                             index=list(range(len(sitelh_ll_list))))
    logging.info(
        "Writing sitelh data to  " + sitelh_csv_path)
    sitelh_df.to_csv(
        sitelh_csv_path, index=False)
    logging.info(
    "Sitelh file is of shape {shape} and stored in {path}".format(shape=sitelh_df.shape, path=sitelh_csv_path))
    curr_msa_stats["sitelh_df"] = sitelh_df
    return sitelh_df


def extract_and_update_RaxML_statistics_from_full_data(curr_msa_stats):
    logging.info("Running RaxML statistics from full data and extracting statistics")
    full_data_path = curr_msa_stats["msa_path"]
    full_data_unique_name = curr_msa_stats["prefix"]
    curr_run_directory = os.path.join(curr_msa_stats["prefix"],"raxml_full_data_results")
    create_or_clean_dir(curr_run_directory)
    run_raxml_on_full_dataset(full_data_path, full_data_unique_name , curr_msa_stats, curr_run_directory)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa_path', action='store', type=str, required=True)
    parser.add_argument('--prefix', action='store', type=str,required=True)
    parser.add_argument('--lasso_training_size', action='store', type=str,required=True)
    args = parser.parse_args()
    curr_msa_stats = {"msa_path": args.msa_path,"prefix": args.prefix,"lasso_training_size": int(args.lasso_training_size)}
    log_file_path = args.prefix+ ".log_file"
    if os.path.exists( log_file_path):
        os.remove(log_file_path)
    logging.basicConfig(filename=args.prefix+ ".log_file", level=LOGGING_LEVEL)
    if os.path.exists(curr_msa_stats["prefix"]):
        logging.error("Folder {} is not empty. Pleese choose a different prefix or empty the current folder.".format(curr_msa_stats["prefix"]))
        return
    os.mkdir(curr_msa_stats["prefix"])
    logging.info("Run parameters: {}".format(curr_msa_stats))
    logging.info("Results will be stored in the folder {}".format(curr_msa_stats["prefix"]))
    try:
        logging.info("Running raxml result on full data:")
        extract_and_update_RaxML_statistics_from_full_data(
            curr_msa_stats)
    except RE_RUN_ON_REDUCED_VERSION:
        extract_and_update_RaxML_statistics_from_full_data(
            curr_msa_stats)
    generate_site_lh_data(curr_msa_stats,curr_msa_stats["lasso_training_size"])
    apply_lasso_on_data_and_update_stats(curr_msa_stats)  # calculating positions_weight



if __name__ == "__main__":
    main()
