import subprocess
import re
from help_functions import *
import os.path
from spr_prune_and_regraft import *
from datetime import datetime



class RE_RUN_ON_REDUCED_VERSION(Exception):
    """Raised when the input value is too large"""
    pass

class GENERAL_RAXML_ERROR(Exception):
    pass

def execute_commnand_and_write_to_log(command_str):
    logging.debug("About to run " + command_str)
    subprocess.run(command_str, shell=True)
    logging.debug("   #Previous command completed")

def check_file_existence(path,name):
    if os.path.exists(path):
        logging.debug("{name} was succesfully created in: {path}".format(name=name,path=path))
    else:
        error_msg="{name} was not generated in: {path}".format(name=name,path=path)
        logging.error(error_msg)
        raise  GENERAL_RAXML_ERROR(error_msg)


def extract_param_from_log(raxml_log_path, param_name):
    with open(raxml_log_path) as raxml_log_file:
        data = raxml_log_file.read()
        if (param_name == "alpha"):
            pattern = r'alpha: ([\d.]+)'
        if (param_name == "ll"):
            pattern = r'Final LogLikelihood: (-[\d.]+)'
        elif (param_name=="rf_dist"):
            pattern='Average relative RF distance in this tree set: ([\d.]+)'
        match = re.search(pattern, data, re.IGNORECASE)
        if match:
                value = float(match.group(1))
        else:
            error_msg= "Param {param_name} not found in file".format(param_name)
            logging.error(error_msg)
            raise GENERAL_RAXML_ERROR(error_msg)
        return value

def extract_mad_file_statistic(mad_log_path):
    pattern  ="MAD=([\d.]+)"
    with open(mad_log_path) as mad_output:
        data = mad_output.read()
        match = re.search(pattern, data, re.IGNORECASE)
    if match:
        value = float(match.group(1))
    else:
        error_msg = "Param  not found in mad file in {}".format(mad_log_path)
        logging.error(error_msg)
        raise GENERAL_RAXML_ERROR(error_msg)
    return value

def calculate_rf_dist(rf_file_path,curr_run_directory):
    rf_prefix = os.path.join(curr_run_directory, "rf")
    rf_command = (
           "{raxml_exe_path} --rfdist --tree {rf_file_path} --prefix {prefix}").format(
        raxml_exe_path =RAXML_NG_COMMAND_PREFIX, rf_file_path=rf_file_path, prefix=rf_prefix)
    execute_commnand_and_write_to_log(rf_command)
    rf_log_file_path = rf_prefix+".raxml.log"
    relative_rf_dist = extract_param_from_log(rf_log_file_path, "rf_dist")
    return relative_rf_dist


def extract_raxml_statistics_from_msa(full_file_path, output_name, msa_stats, curr_run_directory):
    check_validity_prefix = os.path.join(curr_run_directory, output_name + "_CHECK")
    check_validity_command = (
           "{raxml_exe_path} --check --msa {msa_path} --model WAG+G --prefix {prefix}").format(raxml_exe_path =RAXML_NG_COMMAND_PREFIX,
                                                                                               msa_path=full_file_path, prefix=check_validity_prefix)
    reduced_file = check_validity_prefix + ".raxml.reduced.phy"
    execute_commnand_and_write_to_log(check_validity_command)
    if os.path.exists(reduced_file):
        logging.error("Need to re-calculate data on reduced version in " + reduced_file)
        msa_stats["orig_reduced_file_path"] = reduced_file
        raise RE_RUN_ON_REDUCED_VERSION("Input MSA is not valid, re-running on a reduced version")
    parsimony_tree_generation_prefix = os.path.join(curr_run_directory, output_name + "pars")
    parsimony_tree_generation_command = (
           "{raxml_exe_path} --start --msa {msa_path} --model WAG+G --tree pars{{{n_parsimony_trees}}} --prefix {prefix}").format(raxml_exe_path =RAXML_NG_COMMAND_PREFIX,
                                                                                                                                  msa_path=full_file_path, n_parsimony_trees=1, prefix=parsimony_tree_generation_prefix)
    execute_commnand_and_write_to_log(parsimony_tree_generation_command)
    constant_branch_length_parsimony_tree_path = parsimony_tree_generation_prefix + ".raxml.startTree"
    msa_stats["raxml_parsimony_tree_path"]=constant_branch_length_parsimony_tree_path
    check_file_existence(constant_branch_length_parsimony_tree_path, "Parsimony tree")
    parsimony_model_evaluation_prefix = os.path.join(curr_run_directory, output_name + "pars_eval")
    parsimony_model_and_bl_evaluation_command = (
         "{raxml_exe_path} --evaluate --msa {msa_path} --model WAG+G  --tree {parsimony_tree_path} --prefix {prefix}").format(raxml_exe_path =RAXML_NG_COMMAND_PREFIX,
                                                                                                                              msa_path=full_file_path, parsimony_tree_path=constant_branch_length_parsimony_tree_path, prefix=  parsimony_model_evaluation_prefix)
    execute_commnand_and_write_to_log(parsimony_model_and_bl_evaluation_command)
    parsimony_log_path = parsimony_model_evaluation_prefix+".raxml.log"
    check_file_existence(parsimony_log_path, "Parsimony log")
    parsimony_optimized_tree_path=parsimony_model_evaluation_prefix+".raxml.bestTree"
    parsimony_divergence = compute_tree_divergence(parsimony_optimized_tree_path)
    parsimony_tree_alpha = extract_param_from_log(parsimony_log_path, "alpha")
    mad_command = "{mad_exe_path} -t -s {tree_path}".format(mad_exe_path=MAD_COMMAND_PREFIX,tree_path = parsimony_optimized_tree_path)
    execute_commnand_and_write_to_log(mad_command)
    mad_log_path = parsimony_optimized_tree_path+".rooted"
    mad =  extract_mad_file_statistic(mad_log_path)
    msa_stats["mad"]=mad
    msa_stats["alpha"] = parsimony_tree_alpha
    msa_stats["divergence"]=parsimony_divergence


def raxml_extract_sitelh(sitelh_file):
    logging.debug("Extracting sitelh from sitelh_file in {}".format(sitelh_file))
    with open(sitelh_file) as SITELH:
        all_lines = SITELH.readlines()
        sitelh = (re.sub(r'tree1', '', all_lines[1])).strip()
        sitelh_list = sitelh.split(" ")
        # print(sitelh_list)
        sitelh_list_float = [float(ll) for ll in sitelh_list]
        return (sitelh_list_float)



def generate_random_tree_topology(alpha, original_file_path, random_tree_generation_prefix):
    random.seed(datetime.now())
    seed = numpy.random.randint(low=0, high=sys.maxsize)
    random_tree_generation_command = (
            "{raxml_exe_path}  --msa {msa_path} --model WAG+G{{{alpha}}} --start --tree rand{{{n_random_trees}}} --prefix {prefix} --seed {seed}").format(raxml_exe_path =RAXML_NG_COMMAND_PREFIX,
                                                                                                                                            msa_path=original_file_path, alpha=alpha, n_random_trees=1, prefix=random_tree_generation_prefix, seed=seed)
    execute_commnand_and_write_to_log(random_tree_generation_command)
    random_tree_path = random_tree_generation_prefix + ".raxml.startTree"
    check_file_existence(random_tree_path,"random tree")
    return random_tree_path




def raxml_compute_tree_per_site_ll(curr_run_directory, full_data_path, tree_file, ll_on_data_prefix, alpha,opt_brlen=True):
    compute_site_ll_prefix = os.path.join(curr_run_directory, ll_on_data_prefix)
    brlen_command = "--opt-branches off" if not opt_brlen else ""
    compute_site_ll_run_command = (
            "{raxml_exe_path} --sitelh --msa {msa_path} --model WAG+G{{{alpha}}} {brlen_command} --tree {tree_file} --prefix {compute_site_ll_prefix} ").format(raxml_exe_path =RAXML_NG_COMMAND_PREFIX,
                                                                                                                               alpha=alpha, msa_path=full_data_path, tree_file=tree_file,
                                                                                                                               prefix=compute_site_ll_prefix,brlen_command=brlen_command, compute_site_ll_prefix=compute_site_ll_prefix)
    execute_commnand_and_write_to_log( compute_site_ll_run_command)
    sitelh_file = compute_site_ll_prefix + ".raxml.siteLH"
    check_file_existence(sitelh_file,"Sitelh file")
    sitelh_list = raxml_extract_sitelh(sitelh_file)
    return sitelh_list


# def raxml_compute_ll_on_given_data(full_data_path, ll_on_data_prefix, tree_file, msa_stats,
#                                    curr_run_directory, sitelh_lambda_function=None):
#     curr_run_directory = os.path.join(curr_run_directory , ll_on_data_prefix)
#     if os.path.exists(curr_run_directory):
#         delete_dir_content(curr_run_directory)
#     else:
#         os.mkdir(curr_run_directory)
#     logging.debug("RaxML: Evaluating likelihood on : " + full_data_path)
#     alpha = msa_stats["alpha"]
#     logging.debug("Optimizing branch lengths and using existing Gamma shape parameter: alpha={alpha}".format(alpha=alpha))
#     sitelh_list = raxml_compute_tree_per_site_ll(curr_run_directory, full_data_path, tree_file,
#                                                  ll_on_data_prefix, alpha)
#     if sitelh_lambda_function:
#         tree_ll_on_data = sitelh_lambda_function(sitelh_list)
#     else:
#         tree_ll_on_data = sum(sitelh_list)
#     logging.debug("tree ll on data=" + str(tree_ll_on_data))
#     # if os.path.exists(curr_run_directory):
#     #   delete_dir_content(curr_run_directory)
#     return tree_ll_on_data


def raxml_compute_ll_on_given_data(full_data_path, ll_on_data_prefix, tree_file, msa_stats,
                                   curr_run_directory, use_weights=False):
    curr_run_directory = os.path.join(curr_run_directory , ll_on_data_prefix)
    if os.path.exists(curr_run_directory):
        delete_dir_content(curr_run_directory)
    else:
        os.mkdir(curr_run_directory)
    logging.debug("RaxML: Evaluating likelihood on : " + full_data_path)
    alpha = msa_stats["alpha"]
    weights_path_command = "--site-weights {}".format(msa_stats["weights_file_path"]) if use_weights else ""
    logging.debug("Optimizing branch lengths and using existing Gamma shape parameter: alpha={alpha}".format(alpha=alpha))
    prefix = os.path.join(curr_run_directory, ll_on_data_prefix)
    compute_ll_run_command = (
        "{raxml_exe_path} --evaluate --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {weights_path_command} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_COMMAND_PREFIX,
        alpha=alpha, msa_path=full_data_path, tree_file=tree_file,
        prefix=prefix, weights_path_command=weights_path_command)
    execute_commnand_and_write_to_log(compute_ll_run_command)
    raxml_log_file = prefix+".raxml.log"
    tree_ll_on_data = extract_param_from_log(raxml_log_file, "ll")
    logging.debug("tree ll on data=" + str(tree_ll_on_data))
    return tree_ll_on_data

