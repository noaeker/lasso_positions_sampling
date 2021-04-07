
import re
from help_functions import *
import os.path
from spr_prune_and_regraft import *


class RE_RUN_ON_REDUCED_VERSION(Exception):
    """Raised when the input value is too large"""
    pass


class GENERAL_RAXML_ERROR(Exception):
    pass


def execute_commnand_and_write_to_log(command_str):
    logging.debug("About to run " + command_str)
    subprocess.run(command_str, shell=True)
    logging.debug("   #Previous command completed")


def check_file_existence(path, name):
    if os.path.exists(path):
        logging.debug("{name} was succesfully created in: {path}".format(name=name, path=path))
    else:
        error_msg = "{name} was not generated in: {path}".format(name=name, path=path)
        logging.error(error_msg)
        raise GENERAL_RAXML_ERROR(error_msg)

def generate_raxml_command_prefix(cpus=1):
  raxml_parallel_command = " --threads auto{{{N}}} --workers auto ".format(N=cpus)
  return raxml_parallel_command

def extract_param_from_log(raxml_log_path, param_name):
    with open(raxml_log_path) as raxml_log_file:
        data = raxml_log_file.read()
        if (param_name == "alpha"):
            pattern = r'alpha: ([\d.]+)'
        elif (param_name == "search_ll"):
            pattern = r'Final LogLikelihood: (-[\d.]+)'
        elif (param_name == "rf_dist"):
            pattern = 'Average relative RF distance in this tree set: ([\d.]+)'
        elif (param_name == "time"):
            pattern = 'Elapsed time: ([\d.]+)'
        elif (param_name == "ll"):
            pattern = r'Tree #\d+, final logLikelihood: (-[\d.]+)'
            value_strings = re.findall(pattern, data)
            value_floats = [float(ll) for ll in value_strings]
            if len(value_floats) == 1:
                return value_floats[0]
            else:
                logging.info("{} ll values were extracted from log file".format(len(value_floats)))
                return value_floats
        match = re.search(pattern, data, re.IGNORECASE)
        if match:
            value = float(match.group(1))
        else:
            error_msg = "Param {param_name} not found in file".format(param_name=param_name)
            logging.error(error_msg)
            raise GENERAL_RAXML_ERROR(error_msg)
        return value


def extract_mad_file_statistic(mad_log_path):
    pattern = "MAD=([\d.]+)"
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


def raxml_search(curr_run_directory,msa_path, prefix, curr_msa_stats, n_parsimony_trees, n_random_trees,cpus, nodes, weights=None,
                 starting_trees_path=None):
    alpha = curr_msa_stats["alpha"]
    weights_path_command = "--site-weights {}".format(weights) if weights else ""
    if starting_trees_path:
        starting_trees_command = "--tree {}".format(starting_trees_path)
    else:
        starting_trees_command = "--tree pars{{{n_parsimony_trees}}},rand{{{n_random_trees}}}".format(n_parsimony_trees = n_parsimony_trees,
                                                                                                      n_random_trees = n_random_trees)
    search_prefix = os.path.join(curr_run_directory,prefix)
    search_command = (
        "{raxml_exe_path}  {threads_config} --force msa --msa {msa_path} --model WAG+G{{{alpha}}} {starting_trees_command} {weights_path_command} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path= RAXML_NG_EXE,
        threads_config =  generate_raxml_command_prefix(cpus),
        alpha=alpha, msa_path=msa_path, starting_trees_command=starting_trees_command, seed=SEED,
        prefix=search_prefix, weights_path_command=weights_path_command)
    best_tree_path = search_prefix + ".raxml.bestTree"
    raxml_search_starting_tree_path = search_prefix + ".raxml.startTree"
    all_final_trees_path = search_prefix + ".raxml.mlTrees"
    log_file = search_prefix + ".raxml.log"
    if LOCAL_RUN:
        execute_commnand_and_write_to_log(search_command)
    else:
        job_folder = os.path.join(curr_run_directory,"raxml_run_job")
        submit_linux_job("raxml_search", job_folder, search_command, cpus, nodes)
        while not os.path.exists(best_tree_path):
            time.sleep(WAITING_TIME_CSV_UPDATE)
    elapsed_running_time = extract_param_from_log(log_file, 'time')
    best_ll = extract_param_from_log(log_file, 'search_ll')
    return {'best_ll': best_ll, 'best_tree_path': best_tree_path, 'all_final_trees_path': all_final_trees_path,
            'elapsed_running_time': elapsed_running_time, 'starting_trees_path': raxml_search_starting_tree_path}


def raxml_search_pipeline(curr_run_directory,curr_msa_stats, n_parsimony_trees, n_random_trees,standrad_search):
    if standrad_search:
        standard_search_dict = raxml_search(curr_run_directory,curr_msa_stats["local_alignment_path"], "standard", curr_msa_stats,
                                            n_parsimony_trees, n_random_trees,cpus = curr_msa_stats["n_cpus_full"],nodes = curr_msa_stats["n_nodes_full"],
                                            weights=None, starting_trees_path=None)
        results= {'standard_best_ll': standard_search_dict["best_ll"],'standard_best_tree_path': standard_search_dict["best_tree_path"], 'standard_search_elapsed_time':standard_search_dict["elapsed_running_time"],
                  'standard_starting_trees_path': standard_search_dict["starting_trees_path"]}
    else:
        first_phase_dict = raxml_search(curr_run_directory,curr_msa_stats["sampled_alignment_path"], "first_phase", curr_msa_stats, n_parsimony_trees,
                                        n_random_trees,
                                        cpus=curr_msa_stats["n_cpus_Lasso"], nodes=curr_msa_stats["n_nodes_Lasso"],
                                        weights=curr_msa_stats["weights_file_path"], starting_trees_path=curr_msa_stats["standard_starting_trees_path"] if curr_msa_stats["use_raxml_standard_starting_trees"] else None

                                        )
        first_phase_best_true_ll,first_phase_best_tree = raxml_optimize_trees_for_given_msa(
            curr_msa_stats["local_alignment_path"], "first_phase_ll_eval_on_full" ,
           first_phase_dict["best_tree_path"],
            curr_msa_stats, curr_run_directory=curr_run_directory, weights=None)
        results = {
            'lasso_first_phase_best_ll': first_phase_best_true_ll,
            'lasso_first_phase_best_tree': first_phase_dict["best_tree_path"],
            'lasso_first_phase_elapsed_time': first_phase_dict["elapsed_running_time"],

        }
        if curr_msa_stats["do_raxml_lasso_second_phase"]:
            second_phase_dict = raxml_search(curr_run_directory,curr_msa_stats["local_alignment_path"], "second_phase", curr_msa_stats,
                                             n_parsimony_trees, n_random_trees,ncpus = curr_msa_stats["n_cpus_full"],nodes = curr_msa_stats["n_nodes_full"],
                                             weights=None, starting_trees_path=first_phase_dict["all_final_trees_path"])
            results.update(   {'lasso_second_phase_best_ll': second_phase_dict["best_ll"],
                   'lasso_second_phase_best_tree': second_phase_dict["best_tree_path"],
                   'lasso_second_phase_elapsed_time': second_phase_dict["elapsed_running_time"]})

    return results


def calculate_rf_dist(rf_file_path, curr_run_directory):
    rf_prefix = os.path.join(curr_run_directory, "rf")
    rf_command = (
        "{raxml_exe_path} --rfdist --tree {rf_file_path} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE, rf_file_path=rf_file_path, prefix=rf_prefix)
    execute_commnand_and_write_to_log(rf_command)
    rf_log_file_path = rf_prefix + ".raxml.log"
    relative_rf_dist = extract_param_from_log(rf_log_file_path, "rf_dist")
    return relative_rf_dist


def extract_raxml_statistics_from_msa(full_file_path, output_name, msa_stats, curr_run_directory):
    check_validity_prefix = os.path.join(curr_run_directory, output_name + "_CHECK")
    check_validity_command = (
        "{raxml_exe_path} {threads_config} --check --msa {msa_path} --model WAG+G --prefix {prefix} --seed {seed}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_command_prefix(),
        msa_path=full_file_path, prefix=check_validity_prefix, seed=SEED)
    reduced_file = check_validity_prefix + ".raxml.reduced.phy"
    execute_commnand_and_write_to_log(check_validity_command)
    if os.path.exists(reduced_file):
        logging.error("Need to re-calculate data on reduced version in " + reduced_file)
        msa_stats["orig_reduced_file_path"] = reduced_file
        raise RE_RUN_ON_REDUCED_VERSION("Input MSA is not valid, re-running on a reduced version")
    parse_prefix = os.path.join(curr_run_directory, output_name + "_PARSE")
    parse_command = "{raxml_exe_path} {threads_config} --parse --msa {msa_path} --model WAG+G --prefix {prefix} --seed {seed}".format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_command_prefix(),
        msa_path=full_file_path, prefix=parse_prefix, seed=SEED)
    execute_commnand_and_write_to_log(parse_command)
    parsimony_tree_generation_prefix = os.path.join(curr_run_directory, output_name + "pars")
    parsimony_tree_generation_command = (
        "{raxml_exe_path} {threads_config} --start --msa {msa_path} --model WAG+G --tree pars{{{n_parsimony_trees}}} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_command_prefix(),
        msa_path=full_file_path, n_parsimony_trees=1, prefix=parsimony_tree_generation_prefix, seed=SEED)
    execute_commnand_and_write_to_log(parsimony_tree_generation_command)
    constant_branch_length_parsimony_tree_path = parsimony_tree_generation_prefix + ".raxml.startTree"
    msa_stats["raxml_parsimony_tree_path"] = constant_branch_length_parsimony_tree_path
    check_file_existence(constant_branch_length_parsimony_tree_path, "Parsimony tree")
    parsimony_model_evaluation_prefix = os.path.join(curr_run_directory, output_name + "pars_eval")
    parsimony_model_and_bl_evaluation_command = (
        "{raxml_exe_path} {threads_config} --evaluate --msa {msa_path} --model WAG+G  --tree {parsimony_tree_path} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_command_prefix(),
        msa_path=full_file_path, parsimony_tree_path=constant_branch_length_parsimony_tree_path, seed=SEED,
        prefix=parsimony_model_evaluation_prefix)
    execute_commnand_and_write_to_log(parsimony_model_and_bl_evaluation_command)
    parsimony_log_path = parsimony_model_evaluation_prefix + ".raxml.log"
    check_file_existence(parsimony_log_path, "Parsimony log")
    parsimony_optimized_tree_path = parsimony_model_evaluation_prefix + ".raxml.bestTree"
    msa_stats["parsimony_optimized_tree_path"] = parsimony_optimized_tree_path
    parsimony_divergence = compute_tree_divergence(parsimony_optimized_tree_path)
    parsimony_tree_alpha = extract_param_from_log(parsimony_log_path, "alpha")
    mad_command = "{mad_exe_path} -t -s {tree_path}".format(mad_exe_path=MAD_COMMAND_PREFIX,
                                                            tree_path=parsimony_optimized_tree_path)
    execute_commnand_and_write_to_log(mad_command)
    mad_log_path = parsimony_optimized_tree_path + ".rooted"
    mad = extract_mad_file_statistic(mad_log_path)
    msa_stats["mad"] = mad
    msa_stats["alpha"] = parsimony_tree_alpha
    msa_stats["divergence"] = parsimony_divergence


def raxml_extract_sitelh(sitelh_file):
    logging.debug("Extracting sitelh from sitelh_file in {}".format(sitelh_file))
    with open(sitelh_file) as SITELH:
        sitelh_data = SITELH.read()
        pattern = r'tree\d+\s+([-\d.\s]+)'
        sitelh_strings = re.findall(pattern, sitelh_data)
        sitelh_lists = [sitelh_string.split(" ") for sitelh_string in sitelh_strings]
        sitelh_lists_floats = []
        for sitelh_list in sitelh_lists:
            sitelh_lists_floats.append([float(ll) for ll in sitelh_list if len(ll) > 0])
        return (sitelh_lists_floats)


def generate_n_random_tree_topology_constant_brlen(n, alpha, original_file_path, random_tree_generation_prefix,curr_msa_stats, seed):
    random_tree_generation_command = (
        "{raxml_exe_path} {threads_config} --force msa  --msa {msa_path} --model WAG+G{{{alpha}}} --start --tree rand{{{n}}} --prefix {prefix} --opt-branches off --seed {seed} ").format(
        n=n, raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_command_prefix(cpus= curr_msa_stats["n_cpus_training"]),
        msa_path=original_file_path, alpha=alpha, prefix=random_tree_generation_prefix, seed=seed)
    random_tree_path = random_tree_generation_prefix + ".raxml.startTree"
    log_file = random_tree_generation_prefix + ".raxml.log"
    if LOCAL_RUN:
        execute_commnand_and_write_to_log(random_tree_generation_command)
    else:
        job_folder = os.path.join(random_tree_generation_prefix, "raxml_random_tree_generation_job")
        submit_linux_job("rand_top", job_folder, random_tree_generation_command, curr_msa_stats["n_cpus_training"],
                         curr_msa_stats["n_nodes_training"])
        while not os.path.exists(random_tree_path):
            time.sleep(WAITING_TIME_CSV_UPDATE)
    check_file_existence(random_tree_path, "random tree")
    elapsed_running_time = extract_param_from_log(log_file, 'time')
    return random_tree_path,elapsed_running_time


def raxml_compute_tree_per_site_ll(curr_run_directory, full_data_path, tree_file, ll_on_data_prefix, alpha, curr_msa_stats,
                                   opt_brlen=True):
    compute_site_ll_prefix = os.path.join(curr_run_directory, ll_on_data_prefix)
    brlen_command = "--opt-branches off" if not opt_brlen else ""
    compute_site_ll_run_command = (
        "{raxml_exe_path} {threads_config} --force msa --sitelh --msa {msa_path} --model WAG+G{{{alpha}}} {brlen_command} --tree {tree_file} --seed {seed} --prefix {compute_site_ll_prefix} ").format(
        raxml_exe_path=RAXML_NG_EXE, threads_config=generate_raxml_command_prefix(cpus = curr_msa_stats["n_cpus_training"]),
        alpha=alpha, msa_path=full_data_path, tree_file=tree_file, seed=SEED,
        prefix=compute_site_ll_prefix, brlen_command=brlen_command, compute_site_ll_prefix=compute_site_ll_prefix)
    sitelh_file = compute_site_ll_prefix + ".raxml.siteLH"
    log_file = compute_site_ll_prefix + ".raxml.log"
    if LOCAL_RUN:
        execute_commnand_and_write_to_log(compute_site_ll_run_command)
    else:
        job_folder = os.path.join(curr_run_directory, "raxml_ll_eval_job_for_training")
        submit_linux_job("training_opt", job_folder, compute_site_ll_run_command,curr_msa_stats["n_cpus_training"], curr_msa_stats["n_nodes_training"])
        while not os.path.exists(sitelh_file):
            time.sleep(WAITING_TIME_CSV_UPDATE)
    time.sleep(WAITING_TIME_CSV_UPDATE)
    check_file_existence(sitelh_file, "Sitelh file")
    sitelh_list = raxml_extract_sitelh(sitelh_file)
    elapsed_running_time = extract_param_from_log(log_file, 'time')
    return sitelh_list,elapsed_running_time



def raxml_optimize_trees_for_given_msa(full_data_path, ll_on_data_prefix, tree_file, msa_stats,
                                       curr_run_directory, opt_brlen=True, weights=None, return_trees_file=False):
    curr_run_directory = os.path.join(curr_run_directory, ll_on_data_prefix)
    if os.path.exists(curr_run_directory):
        delete_dir_content(curr_run_directory)
    else:
        os.mkdir(curr_run_directory)
    logging.debug("RaxML: Evaluating likelihood on : " + full_data_path)
    alpha = msa_stats["alpha"]
    weights_path_command = "--site-weights {}".format(weights) if weights else ""
    logging.debug(
        "Optimizing branch lengths and using existing Gamma shape parameter: alpha={alpha}".format(alpha=alpha))
    prefix = os.path.join(curr_run_directory, ll_on_data_prefix)
    brlen_command = "--opt-branches off" if not opt_brlen else ""
    compute_ll_run_command = (
        "{raxml_exe_path} {threads_config} --force msa --evaluate --msa {msa_path} --model WAG+G{{{alpha}}} {brlen_command} --tree {tree_file} {weights_path_command} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_command_prefix(msa_stats["n_cpus_full"]),
        alpha=alpha, msa_path=full_data_path, tree_file=tree_file, seed=SEED,
        prefix=prefix, weights_path_command=weights_path_command, brlen_command=brlen_command)
    optimized_trees_path = prefix + ".raxml.mlTrees"
    best_tree_path = prefix + ".raxml.bestTree"
    if LOCAL_RUN:
        execute_commnand_and_write_to_log( compute_ll_run_command)
    else:
        job_folder = os.path.join(curr_run_directory, "raxml_optimize_test_trees_job")
        submit_linux_job("test_opt", job_folder, compute_ll_run_command,msa_stats["n_cpus_training"], msa_stats["n_nodes_training"])
        while not os.path.exists(best_tree_path):
            time.sleep(WAITING_TIME_CSV_UPDATE)
    time.sleep(WAITING_TIME_CSV_UPDATE)
    raxml_log_file = prefix + ".raxml.log"
    trees_ll_on_data = extract_param_from_log(raxml_log_file, "ll")
    optimized_trees_final_path = optimized_trees_path if os.path.exists(optimized_trees_path) else best_tree_path
    tree_objects = generate_multiple_tree_object_from_newick( optimized_trees_final_path)
    if return_trees_file:
        return optimized_trees_path
    return trees_ll_on_data, tree_objects
