
import logging
import random


LOCAL_RUN = False
SEED = 1

N_PARSIMONY_RAXML_SEARCH= 1
N_RANDOM_RAXML_SEARCH = 1
CPUS_PER_NODE = 1
CPUS_PER_NODE_LASSO = 1
N_NODES = 1
N_NODES_LASSO=1
MPI_PROC_PER_NODE = 1


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
LOGGING_LEVEL = logging.DEBUG
GENERATE_LASSO_DESCRIPTIVE = True
RANDOM_TREES_TRAINING_SIZE = 400
TRAINING_SIZE_OPTIONS = [100,200,400,800,1600,3200]
BRLEN_GENERATORS = {'exponential':sample_exp,'uniform': sample_uniform,'optimized': None}

WAITING_TIME_UPDATE = 60 #86400
N_JOBS = 1
RANDOM_TREES_TEST_SIZE = 2
DELETE_SPR_FILES = True
EPSILON = 0.1
MSA_EXTRACTION_METHOD = "CSV"  # MSA_EXTRACTION_METHOD = "FOLDER"

USE_INTEGER_WEIGHTS = LOCAL_RUN
INTEGER_CONST = 1000 if USE_INTEGER_WEIGHTS else 1
CURR_RUN_PREFIX = "lasso_test_cv_search_test1"
CURR_JOBS_PREFIX = "lasso_test_cv_search_test1"
LASSO_BASELINE = "no_baseline"#"new_test" #"test_unbiassed_lasso"#"test_new"#"test_lasso_random" #"raxml_search_test"
TRAINING_BASELINE = "new_test" #"opt_new_tests_30"#"test_alpha"
TEST_SET_BASELINE = "new_test"
MSA_BASELINE = "new_test"
FULL_DATA_BASELINE = "new_test2_raxml" #"opt_new_tests_30"#"test_unbiassed_lasso"#"raxml_search_test_standard"#"spr_baseline"

DILUTE_AMOUNT = 100
DILUTE_MUL = 10

MAX_N_SEQ = 200
MIN_N_SEQ = 5
N_RANDOM_STARTING_TREES = 1
#PARSIMONY_STARTING_TREE = False #1/0
N_MSAS = 1
FIRST_MSA_IND = 0

OUTPUT_CSV_NAME = "spr_raxml"

IGNORE_COLS_IN_CSV = ["alignment_data","MSA_original_alignment_data", "lasso_coeffs", "lasso_chosen_weights", "lasso_chosen_locis","lasso_predict_func"]
MODULE_LOAD_STR = "module load gcc/gcc-8.2.0; module load python/python-anaconda3.6.5-orenavr2; module load intel/parallel_studio_xe_2020.4.omnipath;"

#module load mpi/openmpi-x86_64

if not LOCAL_RUN:
    # PATH CONFIGURATION
    RAXML_NG_EXE = "/groups/pupko/noaeker/raxml-ng-float-mpi/raxml-ng --extra thread-pin "
    MAD_COMMAND_PREFIX = "/groups/pupko/noaeker/mad"
    RESULTS_FOLDER = "/groups/pupko/noaeker/lasso_positions_sampling_results"
    MSAs_FOLDER = "/groups/pupko/noaeker/data/ABC_DR"
    MSAs_CSV_PATH = "/groups/pupko/noaeker/data/sampled_datasets.csv"
    ALTERNATIVER_FILES_FOLDER = "/groups/pupko/noaeker/example"
    MAIN_CODE_PATH = "/groups/pupko/noaeker/lasso_positions_sampling/parallel_code/MSA_positions_sampling.py"
elif LOCAL_RUN:
    IQTREE_PATH = "/Users/noa/Programs/iqtree-1.6.12-MacOSX/iqtree"
    RAXML_NG_EXE = "/Users/noa/Programs/Raxml/raxml-ng  "
    MAD_COMMAND_PREFIX = "/Users/noa/Programs/mad.osx"
    RESULTS_FOLDER= "/Users/noa/Workspace/lasso_positions_sampling_results"
    MSAs_FOLDER = "/Users/noa/Workspace/data/ABC_DR"
    MSAs_CSV_PATH = "/Users/noa/Workspace/data/sampled_datasets.csv"
    ALTERNATIVER_FILES_FOLDER= "/Users/noa/Workspace/data/LARGE_FILES"
    MAIN_CODE_PATH = "/Users/noa/Workspace/lasso_positions_sampling/parallel_code/MSA_positions_sampling.py"



