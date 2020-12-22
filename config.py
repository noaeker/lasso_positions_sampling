import numpy as np
import logging
import os

LINUX_LOCAL = True#True

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
LOGGING_LEVEL = logging.DEBUG
GENERATE_LASSO_DESCRIPTIVE = True
RANDOM_TREES_TRAINING_SIZE = 1000
WAITING_TIME_CSV_UPDATE = 10#86400
NJOBS =50
RANDOM_TREES_TEST_SIZE = 100
DELETE_SPR_FILES = True
DELETE_CSV_IF_EXISTS = True
EPSILON = 0.1
USE_BACKUP_CSV_FILES_IF_EXISTS = False
UPDATE_BACKUP_IF_EXISTS = False
MSA_EXTRACTION_METHOD = "FOLDER"  # MSA_EXTRACTION_METHOD = "FOLDER"

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



if not LINUX_LOCAL:
    # PATH CONFIGURATION
    ALIGNMENT_CSV_PATH ="/groups/pupko/noaeker/positions_sampling/ABC_DR/sampled_datasets.csv"
    BASELINE_FOLDER = "/groups/pupko/noaeker/"
    PROJECT_FOLDER = BASELINE_FOLDER +"positions_sampling/"
    RAXML_SOURCE_CODE_PATH = BASELINE_FOLDER + "standard-RAxML"
    RAXML_NG_EXECUTABLE_PATH = os.path.join(BASELINE_FOLDER , "raxml-ng/raxml-ng --threads 1")
    ALIGNMENTS_FOLDER_PATH = BASELINE_FOLDER + "Alignment_files"
    ALL_RUNS_BACKUP_FOLDER = PROJECT_FOLDER + ("backup_files")
    MSA_FILES_PREFIX_TO_REMOVE = "/groups/pupko/noaeker"

elif LINUX_LOCAL:
    BASELINE_FOLDER ="/home/noa/"
    RAXML_NG_EXECUTABLE_PATH = BASELINE_FOLDER +"raxml-ng/bin/raxml-ng"
    ALIGNMENTS_FOLDER_PATH = BASELINE_FOLDER + "positions_sampling/Alignment_files"
    PROJECT_FOLDER = BASELINE_FOLDER + "positions_sampling"
    ALL_RUNS_BACKUP_FOLDER = PROJECT_FOLDER + ("backup_files")
    MSA_FILES_PREFIX_TO_REMOVE = BASELINE_FOLDER
    MSA_CODE_LOCATION = "/mnt/c/Users/n1234/PycharmProjects/Current_column_sampling/MSA_positions_sampling.py"



