
job_line = 'module load gcc/gcc-8.2.0; module load python/python-anaconda3.6.5-orenavr2!@#python; python /groups/pupko/noaeker/positions_sampling/current_code/MSA_positions_sampling.py --run_prefix five_taxa --jobs_prefix five_taxa --n_MSAs 50 --n_jobs 50 --first_msa_ind 0 --n_random_starting_trees 1000 --max_n_seq 5\tfive_taxa'

cmds_path  = "five_taxa.cmds"

with open(cmds_path,'w') as cmds:
    cmds.write(job_line)
# if not LOCAL_RUN:
#     os.system(f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {cmds_path} {job_log_path}')