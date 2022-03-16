from config import *
import argparse


def main_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--logging_level', type=str, default=LOGGING_LEVEL)
    parser.add_argument('--run_prefix', action='store', type=str, default=CURR_RUN_PREFIX)  # simply a name
    parser.add_argument('--jobs_prefix', action='store', type=str, default=CURR_JOBS_PREFIX)  # simply a name
    parser.add_argument('--n_MSAs', action='store', type=int, default=N_MSAS)  # number of MSAs to run on
    parser.add_argument('--n_jobs', action='store', type=int, default=N_JOBS)  # Number of jobs
    parser.add_argument('--first_msa_ind', action='store', type=int, default=0)
    parser.add_argument('--n_random_starting_trees', action='store', type=int, default=N_RANDOM_STARTING_TREES)
    parser.add_argument('--random_trees_training_size', action='store', type=str, default=RANDOM_TREES_TRAINING_SIZE)
    parser.add_argument('--exp_brlen', action='store_true', default=True)  # use exponential branch lengths
    parser.add_argument('--uni_brlen', action='store_true')
    parser.add_argument('--opt_brlen', action='store_true')
    parser.add_argument('--const_brlen', action='store_true')
    parser.add_argument('--random_trees_test_size', action='store', type=int,
                        default=RANDOM_TREES_TEST_SIZE)  # size of test set for the Lasso approximation
    parser.add_argument('--max_n_seq', action='store', type=str,
                        default=MAX_N_SEQ)  # maximal number of sequences (used to trim MSA)
    parser.add_argument('--min_n_seq', action='store', type=int,
                        default=MIN_N_SEQ)  # minimal number of sequences (to filter out small MSAs)
    parser.add_argument('--max_n_loci', type=str,
                        default=MAX_N_LOCI)  # maximal number of sites in the MSA (used to trim MSA)
    parser.add_argument('--min_n_loci', type=int,
                        default=MIN_N_LOCI)  # minimal number of sites in the MSA (to filter out small MSAs)
    parser.add_argument('--only_evaluate_lasso', action='store_true',
                        default=True)  # only do Lasso analysis without SPR search
    parser.add_argument('--training_set_baseline_run_prefix', action='store', type=str,
                        default=TRAINING_BASELINE)  # for using existing training set
    parser.add_argument('--lasso_baseline_run_prefix', action='store', type=str,
                        default=LASSO_BASELINE)  # for using existing lasso model
    parser.add_argument('--msa_baseline_run_prefix', action='store', type=str,
                        default=MSA_BASELINE)  # for using existing msa data
    parser.add_argument('--spr_baseline_run_prefix', action='store', type=str,
                        default=FULL_DATA_BASELINE)  # for using existing SPR data
    parser.add_argument('--test_set_baseline_run_prefix', action='store', type=str,
                        default=TEST_SET_BASELINE)  # for using existing test data
    parser.add_argument('--RAxML_baseline_run_prefix', action='store', type=str,
                        default=FULL_DATA_BASELINE)  # for using existing RAXML search data
    parser.add_argument('--alternative_training_prefix', action='store', type=str,
                        default=ALTERNATIVE_TRAINING_BASELINE)
    parser.add_argument('--evo_model', type=str, default=EVO_MODEL)  # Evolutionary model (WAG/JTT/GTR...)
    parser.add_argument('--compare_lasso_to_naive',
                        action='store_true')  # Compare Lasso to random sampling and high evolutionary rate sampling
    parser.add_argument('--compare_loci_gene_distribution', action='store_true')  # , default= True
    parser.add_argument('--do_partitioned_lasso_analysis',
                        action='store_true')  # use model files for partitioned analysis
    parser.add_argument('--n_raxml_parsimony_trees', action='store', type=int, default=N_PARSIMONY_RAXML_SEARCH)
    parser.add_argument('--n_raxml_random_trees', action='store', type=int, default=N_RANDOM_RAXML_SEARCH)
    parser.add_argument('--use_raxml_standard_starting_trees', action='store_true', default=True)
    parser.add_argument('--use_raxml_search', action='store_true')  # Do RAxML search instead of a naive SPR search
    parser.add_argument('--queue', type=str, default="pupkolab")
    parser.add_argument('--do_raxml_lasso_nni_optimization', action='store_true')  # change
    parser.add_argument('--alternative_analysis', action='store_true', default=True)
    parser.add_argument('--n_cpus_per_job', action='store', type=int, default=1)  # for running in a cluster
    parser.add_argument('--n_cpus_full', action='store', type=int, default=CPUS_PER_NODE)  # for running in a cluster
    parser.add_argument('--n_cpus_nni', action='store', type=int, default=CPUS_PER_NODE_NNI)  # for running in a cluster
    parser.add_argument('--n_nodes_full', action='store', type=int, default=N_NODES)  # for running in a cluster
    parser.add_argument('--n_cpus_Lasso', action='store', type=int,
                        default=CPUS_PER_NODE_LASSO)  # for running in a cluster
    parser.add_argument('--n_nodes_Lasso', action='store', type=int, default=N_NODES_LASSO)  # for running in a cluster
    parser.add_argument('--n_cpus_training', action='store', type=int,
                        default=CPUS_PER_NODE)  # for running in a cluster
    parser.add_argument('--n_nodes_training', action='store', type=int, default=N_NODES)  # for running in a cluster
    parser.add_argument('--alternative_files_folder', action='store', type=str, default=INPUT_FILES_FOLDER)
    parser.add_argument('--only_full_search', action='store_true')  # CHANGE
    parser.add_argument('--use_parsimony_training_trees', action='store_true', default=False)
    parser.add_argument('--no_test_set', action='store_true')
    parser.add_argument('--n_partitions', type=int, default=1)  # Not relevant
    parser.add_argument('--lasso_thresholds', action='store', type=str,
                        default=LASSO_THRESHOLDS)  # which sampling threshold to use for Lasso evaluation
    parser.add_argument('--lasso_thresholds_search', action='store', type=str, default=THRESHOLDS_TO_USE_DURING_SEARCH)
    parser.add_argument('--run_raxml_commands_locally', action='store_true')
    parser.add_argument('--random_lasso', action='store_true')
    parser.add_argument('--alternative_training', action='store_true')  # training with non-random trees
    parser.add_argument('--use_glmnet_lasso', action='store_true')
    parser.add_argument('--relaxed_lasso', action='store_true')
    parser.add_argument('--use_spr_parsimony_starting_tree', action='store_true')  # ,
    parser.add_argument('--compute_all_true_ll', action='store_true')
    parser.add_argument('--compute_per_site_ll_values', action='store_true')
    parser.add_argument('--top_ind_to_test_per_phase', action='store', type=str, default=TOP_IND_TO_TEST_PER_PHASE)
    parser.add_argument('--loci_shift', action='store', type=int, default=0)
    parser.add_argument('--rearr_dist', type=int, default=5)  # Rearrangment distance for SPR search
    parser.add_argument('--optimized_neighbours_per_iter', type=int, default=-1)
    parser.add_argument('--greedy_SPR', action='store_true')
    parser.add_argument('--use_spr_neighbours_training', action='store_true')  # training based on previous spr trees
    parser.add_argument('--search_epsilon', default=EPSILON)
    parser.add_argument('--use_modified_final_search', action='store_true')
    parser.add_argument('--start_of_starting_tree_ind', type=int, default=0)
    return parser


def job_parser():
    parser = main_parser()
    parser.add_argument('--job_ind', action='store', type=int)
    parser.add_argument('--curr_job_folder', action='store', type=str)
    return parser
