import pickle
from raxml import *

CHOSEN_JOB_AND_FILE = "job_0/_Selectome_Euteleostomi_ENSGT00390000006373_"
#CHOSEN_JOB_AND_FILE = "job_1/_Selectome_Euteleostomi_ENSGT00680000099708_"
PREFIX = '/Users/noa/Workspace/lasso_positions_sampling_results/test_diff_10'

sampled_file_path = '{prefix}/{chosen_job_and_file}/Lasso_folder/optimized/training_400_random_tree_eval/0_fixed_sampled.phy'.format(chosen_job_and_file=CHOSEN_JOB_AND_FILE,prefix= PREFIX)
lasso_weights_file_path = '{prefix}/{chosen_job_and_file}/Lasso_folder/optimized/training_400_random_tree_eval/0_fixed_lasso_weights.txt'.format(chosen_job_and_file=CHOSEN_JOB_AND_FILE,prefix= PREFIX)
full_msa_path = '{prefix}/{chosen_job_and_file}/raxml_full_data_results_0/0_fixed.phy'.format(chosen_job_and_file=CHOSEN_JOB_AND_FILE,prefix= PREFIX)
optimized_tree_path = '{prefix}/{chosen_job_and_file}/raxml_full_data_results_0_fixed/0pars_eval.raxml.bestTree'.format(chosen_job_and_file=CHOSEN_JOB_AND_FILE,prefix= PREFIX)
#'/Users/noa/Workspace/lasso_positions_sampling_results/test_diff/check_lasso_weights.raxml.bestTree'
curr_results_folder = '{prefix}/check'.format(prefix=PREFIX)
lasso_model_dump_path = '{prefix}/{chosen_job_and_file}/Lasso_folder/optimized/training_400_random_tree_eval/lasso_model.sav'.format(chosen_job_and_file=CHOSEN_JOB_AND_FILE,prefix= PREFIX)
msa_stats_dump_path = '{prefix}/{chosen_job_and_file}/curr_msa_stats.dump'.format(chosen_job_and_file=CHOSEN_JOB_AND_FILE, prefix= PREFIX)

OPT_BRLEN = True
brlen_command = "--opt-branches off" if not OPT_BRLEN else ""

with open(lasso_model_dump_path, 'rb') as model:
    lasso_model = pickle.load(model)

intercept=lasso_model.intercept_
if USE_INTEGER_WEIGHTS:
    weights = [int(lasso_model.coef_[ind] * INTEGER_CONST) for ind in range(len(lasso_model.coef_))]
else:
    weights = lasso_model.coef_
lasso_chosen_locis = [ind for ind in range(len(weights)) if weights[ind] != 0]



create_or_clean_dir(curr_results_folder)


with open(lasso_weights_file_path,'r') as weights_f:
    f= weights_f.read().split(" ")
    lasso_weights= [int(n) for n in f if len(n)>0]

msa_length = len(lasso_model.coef_)

with open(msa_stats_dump_path,'rb') as model:
    curr_msa_stats = pickle.load(model)

alpha = curr_msa_stats["alpha"]
prefix_standard =  os.path.join(curr_results_folder, "standard_weights")
prefix_alternative = os.path.join(curr_results_folder, "alternative_weights")
prefix_alternative_make_sure = os.path.join(curr_results_folder, "alternative_weights_eval")
prefix_constant_weights =  os.path.join(curr_results_folder, "double_weights")
prefix_lasso =  os.path.join(curr_results_folder, "check_lasso_weights")
prefix_standard_on_lasso = os.path.join(curr_results_folder, "check_lasso_tree_on_standard")
prefix_lasso_on_standard = os.path.join(curr_results_folder, "check_standard_tree_on_lasso")
prefix_lasso_on_lasso = os.path.join(curr_results_folder, "check_lasso_tree_on_lasso")


compute_standard_ll_run_command = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {brlen_command} --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=full_msa_path, tree_file=optimized_tree_path,
    prefix=prefix_standard, brlen_command= brlen_command)
subprocess.run(compute_standard_ll_run_command, shell=True)
raxml_log_file_standard = prefix_standard + ".raxml.log"
standard_tree = prefix_standard+ ".raxml.bestTree"
tree_ll_on_data_no_weights = extract_param_from_raxmlNG_log(raxml_log_file_standard, "ll")


lasso_weights_path_command = "--site-weights {}".format(lasso_weights_file_path)
compute_lasso_on_standard_ll_run_command = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {weights_path_command}  --opt-branches off --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=sampled_file_path, tree_file=standard_tree,
    prefix=prefix_lasso_on_standard, weights_path_command=lasso_weights_path_command, brlen_command =brlen_command )
subprocess.run(compute_lasso_on_standard_ll_run_command , shell=True)
raxml_log_file_lasso_standard = prefix_lasso_on_standard + ".raxml.log"
lasso_tree = prefix_lasso_on_standard + ".raxml.bestTree"
ll_standard_using_lasso = extract_param_from_raxmlNG_log(raxml_log_file_lasso_standard, "ll") / INTEGER_CONST + intercept







weights_file_path_alternative = "/Users/noa/Workspace/lasso_positions_sampling_results/test_diff/check/weights.check.alternative"
with open(weights_file_path_alternative, 'w') as f:
   f.write(str(2)+" ")
   for i in range(msa_length-1):
       f.write(str(1) + " ")

weights_path_command_alternative = "--site-weights {}".format(weights_file_path_alternative)

compute_ll_run_command_alternative = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {weights_path_command} {brlen_command} --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=full_msa_path, tree_file=optimized_tree_path, weights_path_command = weights_path_command_alternative,brlen_command=brlen_command,
    prefix=prefix_alternative)
subprocess.run(
compute_ll_run_command_alternative , shell=True)
raxml_log_file_alternative = prefix_alternative + ".raxml.log"
alternative_tree = prefix_alternative + ".raxml.bestTree"
tree_ll_on_data_alternative_weights = extract_param_from_raxmlNG_log(raxml_log_file_alternative, "ll")



compute_ll_run_command_alternative_make_sure = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {weights_path_command} --opt-branches off --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=full_msa_path, tree_file=alternative_tree, weights_path_command = weights_path_command_alternative,brlen_command=brlen_command,
    prefix=prefix_alternative_make_sure)
subprocess.run(
compute_ll_run_command_alternative_make_sure  , shell=True)
raxml_log_file_alternative_make_sure = prefix_alternative_make_sure + ".raxml.log"
tree_ll_on_data_alternative_weights_make_sure = extract_param_from_raxmlNG_log(raxml_log_file_alternative_make_sure, "ll")
alternative_tree_make_sure = prefix_alternative_make_sure+ ".raxml.bestTree"




constant_weights_file_path = "/Users/noa/Workspace/lasso_positions_sampling_results/test_diff/check/weights.check"
with open(constant_weights_file_path, 'w') as f:
   for i in range(msa_length):
       f.write(str(200) + " ")

constant_weights_path_command = "--site-weights {}".format(constant_weights_file_path)


compute_constant_weights_ll_run_command = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {weights_path_command} {brlen_command} --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=full_msa_path, tree_file=optimized_tree_path,
    prefix=prefix_constant_weights, weights_path_command=constant_weights_path_command,brlen_command=brlen_command)
subprocess.run(compute_constant_weights_ll_run_command, shell=True)
raxml_log_file_constant_weights = prefix_constant_weights + ".raxml.log"
half_tree_ll_on_data_constant_weights = extract_param_from_raxmlNG_log(raxml_log_file_constant_weights, "ll") / 200
constant_weight_tree = prefix_constant_weights+ ".raxml.bestTree"


sitelh_list = raxml_compute_tree_per_site_ll(curr_results_folder, full_msa_path, optimized_tree_path, "check_full_sitelh", alpha, opt_brlen=OPT_BRLEN)



sitelh_list_partial = raxml_compute_tree_per_site_ll(curr_results_folder, sampled_file_path, optimized_tree_path, "check_sampled_sitelh", alpha, opt_brlen=OPT_BRLEN)
ll_partial_using_weights = (sum([lasso_weights[i]*sitelh_list_partial[i] for i in range(len(sitelh_list_partial))])/INTEGER_CONST)+intercept
partial_tree = os.path.join(curr_results_folder,"check_sampled_sitelh.raxml.bestTree")


ll_partial_based_full_using_weights = (sum([lasso_weights[i]*sitelh_list[j] for i,j in enumerate(lasso_chosen_locis)])/INTEGER_CONST)+intercept


lasso_weights_path_command = "--site-weights {}".format(lasso_weights_file_path)
compute_standard_ll_run_command = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {weights_path_command} {brlen_command} --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=sampled_file_path, tree_file=optimized_tree_path,
    prefix=prefix_lasso, weights_path_command=lasso_weights_path_command, brlen_command =brlen_command )
subprocess.run(compute_standard_ll_run_command, shell=True)
raxml_log_file_lasso = prefix_lasso + ".raxml.log"
lasso_tree = prefix_lasso + ".raxml.bestTree"
ll_based_on_lasso = extract_param_from_raxmlNG_log(raxml_log_file_lasso, "ll") / INTEGER_CONST




compute_lasso_ll_run_command_on_lasso = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {weights_path_command} --opt-branches off --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=sampled_file_path, tree_file=lasso_tree,
    prefix=prefix_lasso_on_lasso  ,weights_path_command=lasso_weights_path_command, brlen_command= brlen_command)
subprocess.run(compute_lasso_ll_run_command_on_lasso, shell=True)
raxml_log_file_lasso_on_lasso =prefix_lasso_on_lasso   + ".raxml.log"
ll_lasso_on_lasso = extract_param_from_raxmlNG_log(raxml_log_file_lasso_on_lasso, "ll") / INTEGER_CONST + intercept


compute_standard_ll_run_command_on_lasso = (
    "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} --opt-branches off --prefix {prefix}").format(
    raxml_exe_path=RAXML_NG_EXE,
    alpha=alpha, msa_path=full_msa_path, tree_file=lasso_tree,
    prefix=prefix_standard_on_lasso , brlen_command= brlen_command)
subprocess.run(compute_standard_ll_run_command_on_lasso, shell=True)
raxml_log_file_standard_on_lasso = prefix_standard_on_lasso  + ".raxml.log"
true_ll_lasso_tree = extract_param_from_raxmlNG_log(raxml_log_file_standard_on_lasso, "ll")



print("original tree:")
get_tree_string(optimized_tree_path)
print("standard ll {} and with constant weights {} and estimated by lasso sampled data and weight {}".format(tree_ll_on_data_no_weights, half_tree_ll_on_data_constant_weights,ll_standard_using_lasso ))
print("standard tree:")
get_tree_string(standard_tree)
print("using constant weights:")
get_tree_string(constant_weight_tree)
print("tree_ll_on_data_alternative= {} make sure using eval = {}".format(tree_ll_on_data_alternative_weights,tree_ll_on_data_alternative_weights_make_sure))
print("alternative tree:")
get_tree_string(alternative_tree)
print("alternative tree make sure:")
get_tree_string(alternative_tree_make_sure)
diff = tree_ll_on_data_alternative_weights-tree_ll_on_data_no_weights
ll_first = sitelh_list[0]
print("diff= {} ll_first = {}".format(diff,ll_first))

print("full with sitleh = ",sum(sitelh_list))
print("Lasso intercept= {} Lasso_likelihood = {}".format(intercept, ll_based_on_lasso))
tree_ll_on_data_lasso = ll_based_on_lasso + intercept
print("overall tree ll on data using Lasso based on weights= {} true likelihood of the obtained tree is {} lasso eval on it is {}".format(tree_ll_on_data_lasso, true_ll_lasso_tree,ll_lasso_on_lasso))
print("Lasso tree:")
get_tree_string(lasso_tree)
print("overall tree ll on data using Lasso by hand, based on partial= {}".format(ll_partial_using_weights))
print("Tree based on the sampled MSA:")
get_tree_string(partial_tree)
print("overall tree ll on data using Lasso by hand, based on full= {}".format(ll_partial_based_full_using_weights))

