import numpy as np
import logging
import os

LOCAL_RUN = True#True

########### GENERAL RUNNING CONFIGURATIONS #################

CSV_COLUMNS = ["job_id","file_ind","dataset_id","alignment_type",
"n_seq","n_seq_before_reduction","original_n_seq","n_loci","original_n_loci",
"raxml_parsimony_tree_path","raxml_parsimony_ll","full_data_likelihood_raxml",
"full_data_running_time_raxml","SPR_chosen_starting_tree_path","SPR_starting_tree_true_ll","full_data_likelihood_SPR","full_data_SPR_moves",
"sample_pct","sample_method","SPR_starting_tree_ll","sampled_data_final_likelihood_SPR",
"sampled_data_final_likelihood_for_each_part_SPR","sampled_data_total_SPR_moves","Using_sampled_data_SPR_moves", "Using_full_data_SPR_moves",
"sampled_data_SPR_moves_for_each_part",
"overall_SPR_neighbours_per_iteration","true_vs_sampled_ll_per_iteration_for_each_part","pearson_sampled_vs._full_MSA_correlation_per_iteration",
"random_trees_sample_size","max_number_of_msa_sequences","starting_tree_type"
               ]

# "sampled data final likelihood for each part raxml","sampled data final likelihood raxml","sampled data total running time raxml","sampled data running time for each part raxml"



#Default values
LOGGING_LEVEL = logging.INFO
GENERATE_LASSO_DESCRIPTIVE = True
RANDOM_TREES_TRAINING_SIZE = 1000
WAITING_TIME_CSV_UPDATE = 10#86400
N_JOBS = 50
RANDOM_TREES_TEST_SIZE = 100
DELETE_SPR_FILES = True
DELETE_CSV_IF_EXISTS = True
EPSILON = 0.1
USE_BACKUP_CSV_FILES_IF_EXISTS = False
UPDATE_BACKUP_IF_EXISTS = False
MSA_EXTRACTION_METHOD = "CSV"  # MSA_EXTRACTION_METHOD = "FOLDER"

CURR_RUN_PREFIX= "test"
CURR_JOBS_PREFIX = "test_job"


MAX_N_SEQ = 20
MIN_N_SEQ = 20
N_RANDOM_STARTING_TREES = 5
PARSIMONY_STARTING_TREE = False #1/0
N_MSAS = 50
FIRST_MSA_IND = 0

DETAILED_SPR_LOG = False
OUTPUT_CSV_NAME = "spr_raxml"

IGNORE_COLS_IN_CSV = ["alignment_data", "lasso_coeffs", "lasso_chosen_weights", "lasso_chosen_locis","lasso_predict_func"]



if not LOCAL_RUN:
    # PATH CONFIGURATION
    RAXML_NG_EXECUTABLE_PATH = "/groups/pupko/noaeker/raxml-ng/raxml-ng"
    RESULTS_FOLDER = "/groups/pupko/noaeker/lasso_positions_sampling_results"
    MSAs_FOLDER = "/groups/pupko/noaeker/data/ABC_DR"
    MSAs_CSV_PATH = "/groups/pupko/noaeker/data/sampled_datasets.csv"
elif LOCAL_RUN:
    RAXML_NG_EXECUTABLE_PATH = "/Users/noa/Programs/Raxml/raxml-ng"
    RESULTS_FOLDER= "/Users/noa/Workspace/lasso_positions_sampling_results"
    MSAs_FOLDER = "/Users/noa/Workspace/data/ABC_DR"
    MSAs_CSV_PATH = "/Users/noa/Workspace/data/sampled_datasets.csv"

MSAs_CSV_PATH = "/groups/pupko/noaeker/sampled_datasets.csv"
MSAs_FILES_PREFIX_TO_REMOVE = ""
MAIN_CODE_PATH = "MSA_positions_sampling.py"

