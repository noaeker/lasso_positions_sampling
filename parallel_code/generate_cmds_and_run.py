from parallel_code.help_functions import *
args = ps_parser()
job_line = f'module load gcc/gcc-8.2.0; module load python/python-anaconda3.6.5-orenavr2!@#python; python /groups/pupko/noaeker/positions_sampling/current_code/MSA_positions_sampling.py --run_prefix {args.run_prefix} --jobs_prefix {args.jobs_prefix} --n_MSAs {args.n_MSAs} --n_jobs {args.n_jobs} --first_msa_ind {args.first_msa_ind} --n_random_starting_trees {args.n_random_starting_trees} --max_n_seq {args.max_n_seq} --only_evaluate_lasso \t{args.jobs_prefix}'
create_dir_if_not_exists(RESULTS_FOLDER)
curent_run_results_folder = os.path.join(RESULTS_FOLDER, args.run_prefix)
create_or_clean_dir(curent_run_results_folder)
cmds_path =  os.path.join(curent_run_results_folder, args.run_prefix + ".cmds")
job_log_path = os.path.join(curent_run_results_folder, args.run_prefix + "_tmp_log")
with open(cmds_path, 'w') as cmds_f:
    cmds_f.write(job_line)
if not LOCAL_RUN:
    os.system(f'/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {cmds_path} {job_log_path}')