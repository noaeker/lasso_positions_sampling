import numpy as np
import logging
import os

LOCAL_RUN =True #True

########### GENERAL RUNNING CONFIGURATIONS #################

#Default values
LOGGING_LEVEL = logging.DEBUG
GENERATE_LASSO_DESCRIPTIVE = True
RANDOM_TREES_TRAINING_SIZE = 50
WAITING_TIME_CSV_UPDATE = 10#86400
N_JOBS = 1
RANDOM_TREES_TEST_SIZE = 50
DELETE_SPR_FILES = True
DELETE_CSV_IF_EXISTS = True
EPSILON = 0.1
USE_BACKUP_CSV_FILES_IF_EXISTS = False
UPDATE_BACKUP_IF_EXISTS = False
MSA_EXTRACTION_METHOD = "CSV"  # MSA_EXTRACTION_METHOD = "FOLDER"
N_THREADS = 1

CURR_RUN_PREFIX= "test_brlen"
CURR_JOBS_PREFIX = "test_brlen"


MAX_N_SEQ = 5
MIN_N_SEQ = 5
N_RANDOM_STARTING_TREES = 1
PARSIMONY_STARTING_TREE = False #1/0
N_MSAS = 1
FIRST_MSA_IND = 0

DETAILED_SPR_LOG = False
OUTPUT_CSV_NAME = "spr_raxml"

IGNORE_COLS_IN_CSV = ["alignment_data","MSA_original_alignment_data", "lasso_coeffs", "lasso_chosen_weights", "lasso_chosen_locis","lasso_predict_func"]

N_THREADS = 1



if not LOCAL_RUN:
    # PATH CONFIGURATION
    RAXML_NG_COMMAND_PREFIX = "/groups/pupko/noaeker/raxml-ng/raxml-ng --threads {} ".format(N_THREADS)
    MAD_COMMAND_PREFIX = "/groups/pupko/noaeker/mad"
    RESULTS_FOLDER = "/groups/pupko/noaeker/lasso_positions_sampling_results"
    MSAs_FOLDER = "/groups/pupko/noaeker/data/ABC_DR"
    MSAs_CSV_PATH = "/groups/pupko/noaeker/data/sampled_datasets.csv"
elif LOCAL_RUN:
    RAXML_NG_COMMAND_PREFIX = "/Users/noa/Programs/Raxml/raxml-ng --threads {} ".format(N_THREADS)
    MAD_COMMAND_PREFIX = "/Users/noa/Programs/mad.osx"
    RESULTS_FOLDER= "/Users/noa/Workspace/lasso_positions_sampling_results"
    MSAs_FOLDER = "/Users/noa/Workspace/data/ABC_DR"
    MSAs_CSV_PATH = "/Users/noa/Workspace/data/sampled_datasets.csv"

MAIN_CODE_PATH = "MSA_positions_sampling.py"

