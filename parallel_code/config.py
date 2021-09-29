
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
CPUS_PER_NODE_NNI = 2


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
RANDOM_TREES_TRAINING_SIZE = "1000"
BRLEN_GENERATORS = {'exponential':sample_exp,'uniform': sample_uniform,'optimized': None}

WAITING_TIME_UPDATE = 60 #86400
N_JOBS = 1
RANDOM_TREES_TEST_SIZE = 30
DELETE_SPR_FILES = True
EPSILON = 0.1
ALPHA_EPSILON = 0.0001
LASSO_THRESHOLDS  = "0.05"
THRESHOLDS_TO_USE_DURING_SEARCH = "0.05"
MSA_EXTRACTION_METHOD = "CSV"  # MSA_EXTRACTION_METHOD = "FOLDER"

USE_INTEGER_WEIGHTS = LOCAL_RUN
INTEGER_CONST = 10 if USE_INTEGER_WEIGHTS else 1
CURR_RUN_PREFIX = "test_curr_newest"
CURR_JOBS_PREFIX =  "test_curr_newest"

LASSO_BASELINE ="x"
TRAINING_BASELINE =    "x"
TEST_SET_BASELINE = "x"
MSA_BASELINE =  "x"
FULL_DATA_BASELINE = "x"
ALTERNATIVE_TRAINING_BASELINE = "x"



MAX_N_SEQ =  50
MAX_N_LOCI = 50000
MIN_N_SEQ = 6
N_RANDOM_STARTING_TREES = 1
#PARSIMONY_STARTING_TREE = False #1/0
N_MSAS = 1
FIRST_MSA_IND = 1

OUTPUT_CSV_NAME = "spr_raxml"

IGNORE_COLS_IN_CSV = ["alignment_data","MSA_original_alignment_data", "lasso_coeffs", "lasso_chosen_weights", "lasso_chosen_locis","lasso_predict_func","lasso_first_phase_ml_trees_objects","orig_reduced_file_path","alternative_analysis","rate4site_scores",]
MODULE_LOAD_STR = "module load gcc/gcc-8.2.0; module load R/3.6.1; module load python/python-anaconda3.6.5-orenavr2; module load intel/parallel_studio_xe_2020.4.omnipath;"


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
    R_CODE_PATH = "/groups/pupko/noaeker/lasso_positions_sampling/R_code/lasso_glmnet.R"
    RAXML_HPC_EXE = "/groups/pupko/noaeker/standard-RAxML/raxmlHPC"
    RATE4SITE_COMMAND_PREFIX = "/groups/pupko/noaeker/rate4site/rate4site"
elif LOCAL_RUN:
    IQTREE_PATH = "/Users/noa/Programs/iqtree-1.6.12-MacOSX/iqtree"
    RAXML_NG_EXE = "/Users/noa/Programs/Raxml/raxml-ng  "
    RAXML_HPC_EXE = "/Users/noa/Programs/standard-RAxML/raxmlHPC-PTHREADS "
    MAD_COMMAND_PREFIX = "/Users/noa/Programs/mad.osx"
    RATE4SITE_COMMAND_PREFIX = "/Users/noa/Programs/rate4site"
    RESULTS_FOLDER= "/Users/noa/Workspace/lasso_positions_sampling_results"
    MSAs_FOLDER = "/Users/noa/Workspace/data/ABC_DR"
    MSAs_CSV_PATH = "/Users/noa/Workspace/data/sampled_datasets.csv"
    ALTERNATIVER_FILES_FOLDER= "/Users/noa/Workspace/data/supermatrices_edited"
    MAIN_CODE_PATH = "/Users/noa/Workspace/lasso_positions_sampling/parallel_code/MSA_positions_sampling.py"
    R_CODE_PATH = "/Users/noa/Workspace/lasso_positions_sampling/R_code/lasso_glmnet.R"



