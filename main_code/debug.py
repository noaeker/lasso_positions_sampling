import pickle


curr_msa_version_stats_dump_baseline = "/Users/noa/Workspace/lasso_positions_sampling_results/dna_lasso30_test_lasso_par/job_0/_Users_noa_Workspace_data_supermatrices_edited_DNA__supermatrices_PrumD6_PrumD6.fasta/n_seq_30/n_loci_20000/curr_msa_stats.dump"
with open(curr_msa_version_stats_dump_baseline, 'rb') as handle:
            curr_msa_stats = pickle.load(handle)
            print(curr_msa_stats["msa_corrected_model_partition_optimized"])
            print(curr_msa_stats["partition_ind_to_name_optimized"])