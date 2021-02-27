
from raxml import *
import pandas as pd
from matplotlib import pyplot as plt

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
        curr_msa_results.loc[:,"overall_best_topology_first_phase"] = best_naive_spr_tree_newick
    else:
        curr_msa_results.loc[:,"overall_best_topology_first_phase"] = best_lasso_spr_first_phase_tree_newick
    if best_naive_spr_ll > best_lasso_spr_second_phase_ll:
        curr_msa_results.loc[:,"overall_best_topology_second_phase"] = best_naive_spr_tree_newick
    else:
        curr_msa_results.loc[:,"overall_best_topology_second_phase"] = best_lasso_spr_second_phase_tree_newick

    curr_msa_results.loc[:,"rf_naive_vs_overall_best_first_phase"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_naive_vs_best_first_phase(row, curr_run_directory),
        axis=1)
    curr_msa_results.loc[:,"rf_naive_vs_overall_best_second_phase"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_naive_vs_best_second_phase(row, curr_run_directory),
        axis=1)
    curr_msa_results.loc[:,"rf_first_phase_vs_overall_best"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_first_phase_vs_best(row, curr_run_directory),
        axis=1)
    curr_msa_results.loc[:,"rf_second_phase_vs_overall_best"] = curr_msa_results.apply(
        lambda row: rel_rf_dist_second_phase_vs_best(row, curr_run_directory),
        axis=1)
    return curr_msa_results





unified_df_path = ("/Users/noa/Workspace/lasso_positions_sampling_results/sp_c.csv")
output_data= pd.DataFrame(
)
unified_df = pd.read_csv(unified_df_path)
curr_run_directory = "/Users/noa/Workspace/lasso_positions_sampling_results/test_dir"
create_dir_if_not_exists(curr_run_directory)
for dataset_id in unified_df["dataset_id"].unique():
    dataset_id_res = unified_df.loc[unified_df["dataset_id"]==dataset_id]
    enrich_curr_msa_results(dataset_id_res,curr_run_directory)
    if output_data.empty:
        output_data= dataset_id_res
    else:
        output_data=pd.concat([output_data,dataset_id_res])
output_data.to_csv("sp_c.csv")






