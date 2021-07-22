
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
LOGGING_LEVEL = logging.INFO
GENERATE_LASSO_DESCRIPTIVE = True
RANDOM_TREES_TRAINING_SIZE = "800"
BRLEN_GENERATORS = {'exponential':sample_exp,'uniform': sample_uniform,'optimized': None}

WAITING_TIME_UPDATE = 60 #86400
N_JOBS = 1
RANDOM_TREES_TEST_SIZE = 30
DELETE_SPR_FILES = True
EPSILON = 0.1
ALPHA_EPSILON = 0.0001
LASSO_THRESHOLDS  = "0.005_0.01_0.05_0.1_0.15"
THRESHOLDS_TO_USE_DURING_SEARCH = "0.005_0.01_0.05_0.1"
MSA_EXTRACTION_METHOD = "CSV"  # MSA_EXTRACTION_METHOD = "FOLDER"

USE_INTEGER_WEIGHTS = LOCAL_RUN
INTEGER_CONST = 1000 if USE_INTEGER_WEIGHTS else 1
CURR_RUN_PREFIX = "test_lasso_sklearn_7_taxa_800_"#"test_various_thresholds_8"
CURR_JOBS_PREFIX =  "test_lasso_sklearn_7_taxa_800_"#"test_various_thresholds_8"
LASSO_BASELINE = "test_lasso_sklearn_7_taxa_800_"#"test_various_thresholds_6"#"new_test" #"test_unbiassed_lasso"#"test_new"#"test_lasso_random" #"raxml_search_test"
TRAINING_BASELINE =    "test_lasso_glmnet_7_taxa_800"#test_various_thresholds_6"#"test_various_thresholds"#"test_20_for_ppt_400"#"no_baseline"#"new_test" #"opt_new_tests_30"#"test_alpha"
TEST_SET_BASELINE =  "test_lasso_glmnet_7_taxa_800"#"test_various_thresholds_6"#"test_various_thresholds"#"test_20_for_ppt_400"#"raxml_results_30_sample_0.1"
MSA_BASELINE =  "test_lasso_glmnet_7_taxa"#"test_various_thresholds_6"#"test_20_for_ppt_400"#"raxml_results_30_sample_0.1"
FULL_DATA_BASELINE = "test_lasso_glmnet_7_taxa"#"test_10_for_ppt"#"raxml_results_30_sample_0.1"#"new_test2_raxml" #"opt_new_tests_30"#"test_unbiassed_lasso"#"raxml_search_test_standard"#"spr_baseline"
ALTERNATIVE_TRAINING_BASELINE = "no_baseline"

DILUTE_AMOUNT = 15
DILUTE_MUL = 10



MAX_N_SEQ =  7
MAX_N_LOCI = 100000000
MIN_N_SEQ = 7
N_RANDOM_STARTING_TREES = 1
#PARSIMONY_STARTING_TREE = False #1/0
N_MSAS = 1
FIRST_MSA_IND = 1

OUTPUT_CSV_NAME = "spr_raxml"

IGNORE_COLS_IN_CSV = ["alignment_data","MSA_original_alignment_data", "lasso_coeffs", "lasso_chosen_weights", "lasso_chosen_locis","lasso_predict_func","lasso_first_phase_ml_trees_objects","orig_reduced_file_path","alternative_analysis"]
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
elif LOCAL_RUN:
    IQTREE_PATH = "/Users/noa/Programs/iqtree-1.6.12-MacOSX/iqtree"
    RAXML_NG_EXE = "/Users/noa/Programs/Raxml/raxml-ng  "
    MAD_COMMAND_PREFIX = "/Users/noa/Programs/mad.osx"
    RESULTS_FOLDER= "/Users/noa/Workspace/lasso_positions_sampling_results"
    MSAs_FOLDER = "/Users/noa/Workspace/data/ABC_DR"
    MSAs_CSV_PATH = "/Users/noa/Workspace/data/sampled_datasets.csv"
    ALTERNATIVER_FILES_FOLDER= "/Users/noa/Workspace/data/LARGE_FILES/Borowiek_et_al_2015"
    MAIN_CODE_PATH = "/Users/noa/Workspace/lasso_positions_sampling/parallel_code/MSA_positions_sampling.py"
    R_CODE_PATH = "/Users/noa/Workspace/lasso_positions_sampling/R_code/lasso_glmnet.R"



