
import logging
import random


LOCAL_RUN = True #True
SEED = 1

RAxML_SEARCH = True
N_PARSIMONY_RAXML_SEARCH= 1
N_RANDOM_RAXML_SEARCH = 1
RAXML_USE_STANDARD_STARTING_TREES = True
DO_RAXML_SECOND_PHASE = False

########### GENERAL RUNNING CONFIGURATIONS #################

def sample_uniform(size,start_seed):
    res=[]
    seed = start_seed
    for i in range(size):
        random.seed(seed)
        res.append(random.random())
        seed=seed+1
    return res

def sample_exp(size,start_seed):
    res = []
    seed = start_seed
    for i in range(size):
        random.seed(seed)
        res.append(random.expovariate(lambd=10))
        seed = seed + 1
    return res



#Default values
LOGGING_LEVEL = logging.INFO
GENERATE_LASSO_DESCRIPTIVE = True
RANDOM_TREES_TRAINING_SIZE = -1
TRAINING_SIZE_OPTIONS = [400]#[100,200,400,800,1600,3200]
BRLEN_GENERATORS = {'exponential':sample_exp,'uniform': sample_uniform,'optimized': None}


WAITING_TIME_CSV_UPDATE = 10#86400
N_JOBS = 1
RANDOM_TREES_TEST_SIZE = 30
DELETE_SPR_FILES = True
EPSILON = 0.1
MSA_EXTRACTION_METHOD = "CSV"  # MSA_EXTRACTION_METHOD = "FOLDER"
N_THREADS = 1

USE_INTEGER_WEIGHTS = LOCAL_RUN
INTEGER_CONST = 100000 if USE_INTEGER_WEIGHTS else 1
CURR_RUN_PREFIX = "raxml_new"
CURR_JOBS_PREFIX = "job_prefix"
LASSO_BASELINE="run_prefix_raxml_test"

SPR_BASELINE="no_baseline"#"spr_baseline"


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
    RAXML_NG_COMMAND_PREFIX = "/groups/pupko/noaeker/raxml-ng-float/raxml-ng --threads {} --force perf_threads ".format(N_THREADS)
    MAD_COMMAND_PREFIX = "/groups/pupko/noaeker/mad"
    RESULTS_FOLDER = "/groups/pupko/noaeker/lasso_positions_sampling_results"
    MSAs_FOLDER = "/groups/pupko/noaeker/data/ABC_DR"
    MSAs_CSV_PATH = "/groups/pupko/noaeker/data/sampled_datasets.csv"
elif LOCAL_RUN:
    RAXML_NG_COMMAND_PREFIX = "/Users/noa/Programs/Raxml/raxml-ng "
    MAD_COMMAND_PREFIX = "/Users/noa/Programs/mad.osx"
    RESULTS_FOLDER= "/Users/noa/Workspace/lasso_positions_sampling_results"
    MSAs_FOLDER = "/Users/noa/Workspace/data/ABC_DR"
    MSAs_CSV_PATH = "/Users/noa/Workspace/data/sampled_datasets.csv"

MAIN_CODE_PATH = "MSA_positions_sampling.py"

