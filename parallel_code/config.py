import numpy as np
import logging
import os
import help_functions

LOCAL_RUN = False #True

########### GENERAL RUNNING CONFIGURATIONS #################

def sample_uniform(size):
    return np.random.uniform(size=size)
def sample_exp(size):
    return np.random.exponential(scale=0.1, size=size)


#Default values
LOGGING_LEVEL = logging.INFO
GENERATE_LASSO_DESCRIPTIVE = True
RANDOM_TREES_TRAINING_SIZE = -1
TRAINING_SIZE_OPTIONS = [100,200,400,800,1600,3200]
BRLEN_GENERATORS = {'exponential':sample_exp,'uniform': sample_uniform,'optimized': None}


WAITING_TIME_CSV_UPDATE = 10#86400
N_JOBS = 1
RANDOM_TREES_TEST_SIZE = 50
DELETE_SPR_FILES = True
EPSILON = 0.1
MSA_EXTRACTION_METHOD = "CSV"  # MSA_EXTRACTION_METHOD = "FOLDER"
N_THREADS = 1

INTEGER_CONST = 100000
CURR_RUN_PREFIX = "using_baseline_test"
CURR_JOBS_PREFIX = "using_baseline_test"
BASELINE="no_baseline"

MAX_N_SEQ = 5
MIN_N_SEQ = 5
N_RANDOM_STARTING_TREES = 2
PARSIMONY_STARTING_TREE = False #1/0
N_MSAS = 1
FIRST_MSA_IND = 0

DETAILED_SPR_LOG = False
OUTPUT_CSV_NAME = "spr_raxml"

IGNORE_COLS_IN_CSV = ["alignment_data","MSA_original_alignment_data", "lasso_coeffs", "lasso_chosen_weights", "lasso_chosen_locis","lasso_predict_func"]

N_THREADS = 1

USE_PREVIOUS_LASSO_RESULTS = True



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

