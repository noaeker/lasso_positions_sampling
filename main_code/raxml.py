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


def edit_num_locis_in_model_file(model_path, n_loci):
    with open(model_path) as MODEL_FILE:
        text = MODEL_FILE.read()
    new_text = re.sub('noname = 1\-\d+',f'noname = 1-{n_loci}',text)
    with open(model_path,'w') as MODEL_FILE:
        MODEL_FILE.write(new_text)



def execute_commnand_and_write_to_log(command, curr_run_directory="", job_folder_name="", job_name="", log_file_path="",
                                      cpus=-1, nodes=-1, queue="pupkolab", extra_file_path="", run_locally=False):
    if LOCAL_RUN or run_locally:
        logging.debug("*** About to run locally " + command)
        subprocess.run(command, shell=True)
        # logging.info("*** Previous command completed")
    else:
        job_folder = os.path.join(curr_run_directory, job_folder_name)
        submit_linux_job(job_name, job_folder, command, cpus, queue=queue)
        logging.debug(f"*** Waiting for elapsed time in log file {log_file_path}")
        while not (os.path.exists(log_file_path) and (
                os.path.exists(extra_file_path) or extra_file_path == "") and (
                           extract_param_from_raxmlNG_log(log_file_path,
                                                          'time',
                                                          raise_error=False) or extract_param_from_raxmlHPC_log(
                       log_file_path,
                       'time',
                       raise_error=False)) is not None):
            time.sleep(30)
        logging.info("*** current time: {} previous job is completed!!***".format(datetime.now()))


def wait_for_file_existence(path, name):
    if not os.path.exists(path):
        # logging.info("{name} was succesfully created in: {path}".format(name=name, path=path))
        error_msg = "{name} was not generated in: {path}".format(name=name, path=path)
        logging.error(error_msg)
        start_time = time.time()
        while not os.path.exists(path):
            time.sleep(WAITING_TIME_UPDATE)
            logging.info("current time {}: file {} does not exist yet in path {}".format(datetime.now(), name, path))
            time.sleep(WAITING_TIME_UPDATE)
            if time.time() - start_time > 3600*24:
                logging.info("Waiting to much for param {}, breaking".format(name))
                break
        raise GENERAL_RAXML_ERROR(error_msg)


def generate_raxml_ng_command_prefix(cpus=1):
    raxml_parallel_command = " --threads {N} --workers auto ".format(
        N=cpus)  # " --threads auto{{{N}}} --workers auto ".format(N=cpus)
    return raxml_parallel_command


def extract_param_from_raxmlHPC_log(raxml_log_path, param_name, raise_error=True):
    with open(raxml_log_path) as raxml_log_file:
        data = raxml_log_file.read()
        if (param_name == "time"):
            pattern = 'Total execution time: ([\d.]+)'
        match = re.search(pattern, data, re.IGNORECASE)
        if match:
            value = float(match.group(1))
            return value
        else:
            error_msg = "Param {param_name} not found in file".format(param_name=param_name)
            if raise_error:
                raise GENERAL_RAXML_ERROR(error_msg)
            else:
                return None


def extract_param_from_raxmlNG_log(raxml_log_path, param_name, raise_error=True):
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
                # logging.info("{} ll values were extracted from log file".format(len(value_floats)))
                return value_floats
        match = re.search(pattern, data, re.IGNORECASE)
        if match:
            value = float(match.group(1))
            return value
        else:
            error_msg = "Param {param_name} not found in file".format(param_name=param_name)
            if raise_error:
                raise GENERAL_RAXML_ERROR(error_msg)
            else:
                return None




def parse_raxml_partition_file(model_file, orig_n_loci):
    locis_partition = np.zeros(orig_n_loci)
    with open(model_file) as MODEL_FILE:
        partitions = MODEL_FILE.readlines()
    for i,raw_partition in enumerate(partitions):
        site_partitions  = re.findall('\s(\d+)\-(\d+)', raw_partition)
        for site_group in site_partitions:
            start =int(site_group[0])-1
            end = int(site_group[1])
            np.put(locis_partition,ind =np.array(range(start,end)),v  = i+1)
    return locis_partition.astype(int)







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


def raxml_search(curr_run_directory, msa_path, prefix, curr_msa_stats, n_parsimony_trees, n_random_trees, cpus, nodes,
                 weights=None,
                 starting_trees_path=None):
    alpha = curr_msa_stats["alpha"]
    weights_path_command = "--site-weights {}".format(weights) if weights else ""
    if starting_trees_path:
        starting_trees_command = "--tree {}".format(starting_trees_path)
    else:
        starting_trees_command = "--tree pars{{{n_parsimony_trees}}},rand{{{n_random_trees}}}".format(
            n_parsimony_trees=n_parsimony_trees,
            n_random_trees=n_random_trees)
    search_prefix = os.path.join(curr_run_directory, prefix)
    search_command = (
        "{raxml_exe_path}  {threads_config} --force msa --force perf_threads --msa {msa_path} --model WAG+G{{{alpha}}} {starting_trees_command} {weights_path_command} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_ng_command_prefix(cpus),
        alpha=alpha, msa_path=msa_path, starting_trees_command=starting_trees_command, seed=SEED,
        prefix=search_prefix, weights_path_command=weights_path_command)
    best_tree_path = search_prefix + ".raxml.bestTree"
    raxml_search_starting_tree_path = search_prefix + ".raxml.startTree"
    all_final_trees_path = search_prefix + ".raxml.mlTrees" if n_random_trees + n_parsimony_trees > 1 else best_tree_path
    raxml_log_file = search_prefix + ".raxml.log"
    execute_commnand_and_write_to_log(search_command, curr_run_directory, job_folder_name="raxml_search_job",
                                      job_name="raxml_search", log_file_path=raxml_log_file,
                                      cpus=cpus, nodes=nodes, queue=curr_msa_stats["queue"],
                                      run_locally=curr_msa_stats["run_raxml_commands_locally"])
    elapsed_running_time = extract_param_from_raxmlNG_log(raxml_log_file, 'time')
    best_ll = extract_param_from_raxmlNG_log(raxml_log_file, 'search_ll')
    return {'best_ll': best_ll, 'best_tree_path': best_tree_path, 'all_final_trees_path': all_final_trees_path,
            'elapsed_running_time': elapsed_running_time, 'starting_trees_path': raxml_search_starting_tree_path}


def raxml_search_nni(curr_run_directory, msa_path, starting_tree_path, curr_msa_stats, cpus):
    name = "nni_search"
    search_command = (
        "{raxml_exe_path} -T {cpus} -s {msa_path} -f J -m PROTGAMMAWAG -t {starting_tree_path} -p {seed} -w {curr_run_directory} -n {name}").format(
        raxml_exe_path=RAXML_HPC_EXE,
        msa_path=msa_path, starting_tree_path=starting_tree_path, curr_run_directory=curr_run_directory, name=name,
        seed=SEED, cpus=cpus)
    logging.info("search_command is {command}".format(command=search_command))
    best_tree_path = os.path.join(curr_run_directory, "RAxML_fastTree.{name}".format(name=name))
    raxml_log_file = os.path.join(curr_run_directory, "RAxML_info.{name}".format(name=name))
    execute_commnand_and_write_to_log(search_command, curr_run_directory, job_folder_name="raxml_search_job",
                                      job_name="raxml_search", log_file_path=raxml_log_file,
                                      cpus=cpus, nodes=1, queue=curr_msa_stats["queue"],
                                      run_locally=curr_msa_stats["run_raxml_commands_locally"])
    nni_elapsed_running_time = extract_param_from_raxmlHPC_log(raxml_log_file, 'time')
    return {'best_tree_path': best_tree_path,
            'elapsed_running_time': nni_elapsed_running_time}


def raxml_search_pipeline(curr_run_directory, curr_msa_stats, n_parsimony_trees, n_random_trees, standrad_search):
    if standrad_search:
        standard_search_dict = raxml_search(curr_run_directory, curr_msa_stats["local_alignment_path"], "standard",
                                            curr_msa_stats,
                                            n_parsimony_trees, n_random_trees, cpus=curr_msa_stats["n_cpus_full"],
                                            nodes=curr_msa_stats["n_nodes_full"],
                                            weights=None, starting_trees_path=None)

        results = {'standard_best_ll': standard_search_dict["best_ll"],
                   'standard_ml_trees_path': standard_search_dict["all_final_trees_path"],
                   'standard_best_tree_path': standard_search_dict["best_tree_path"],
                   'standard_search_elapsed_time': standard_search_dict["elapsed_running_time"],
                   'standard_starting_trees_path': standard_search_dict["starting_trees_path"]}
    else:
        weights_file_path = curr_msa_stats["weights_file_path"]
        sampled_alignment_path = curr_msa_stats["sampled_alignment_path"]
        logging.info("Performing search using Lasso")

        first_phase_dict = raxml_search(curr_run_directory, sampled_alignment_path, "first_phase",
                                        curr_msa_stats, n_parsimony_trees,
                                        n_random_trees,
                                        cpus=curr_msa_stats["n_cpus_Lasso"], nodes=curr_msa_stats["n_nodes_Lasso"],
                                        weights=weights_file_path,
                                        starting_trees_path=curr_msa_stats["standard_starting_trees_path"] if
                                        curr_msa_stats["use_raxml_standard_starting_trees"] else None

                                        )
        first_phase_best_true_ll, first_phase_best_tree, elapsed_running_time = raxml_optimize_trees_for_given_msa(
            curr_msa_stats["local_alignment_path"], "first_phase_ll_eval_on_full",
            first_phase_dict["best_tree_path"],
            curr_msa_stats, curr_run_directory=curr_run_directory, weights=None, return_trees_file=True)
        results = {
            'lasso_first_phase_best_ll': first_phase_best_true_ll,
            'lasso_first_phase_best_tree': first_phase_dict["best_tree_path"],
            'lasso_first_phase_elapsed_time': first_phase_dict["elapsed_running_time"],
            'lasso_first_phase_ml_trees': first_phase_dict["all_final_trees_path"]

        }
        if curr_msa_stats["do_raxml_lasso_nni_optimization"]:
            logging.info("Performing extra nni optimization")
            nni_dict = raxml_search_nni(curr_run_directory, msa_path=curr_msa_stats["local_alignment_path"],
                                        starting_tree_path=first_phase_best_tree, curr_msa_stats=curr_msa_stats,
                                        cpus=curr_msa_stats["n_cpus_nni"])
            nni_best_ll, nni_best_tree, nni_eval_elapsed = raxml_optimize_trees_for_given_msa(
                curr_msa_stats["local_alignment_path"], "nni_eval",
                nni_dict["best_tree_path"],
                curr_msa_stats, curr_run_directory=curr_run_directory, weights=None, return_trees_file=True)
            results.update({'lasso_nni_best_ll': nni_best_ll,
                            'lasso_nni_best_tree': nni_best_tree,
                            'lasso_nni_elapsed_time': nni_dict["elapsed_running_time"]})

    return results


def calculate_rf_dist(rf_file_path, curr_run_directory, prefix="rf"):
    rf_prefix = os.path.join(curr_run_directory, prefix)
    rf_command = (
        "{raxml_exe_path} --force msa --force perf_threads --rfdist --tree {rf_file_path} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE, rf_file_path=rf_file_path, prefix=rf_prefix)
    execute_commnand_and_write_to_log(rf_command, run_locally=True)
    rf_log_file_path = rf_prefix + ".raxml.log"
    relative_rf_dist = extract_param_from_raxmlNG_log(rf_log_file_path, "rf_dist")
    return relative_rf_dist


def rf_distance(curr_run_directory, tree_object_a, tree_object_b, name):
    rf_folder = os.path.join(curr_run_directory, f"rf_calculations_{name}")
    create_dir_if_not_exists(rf_folder)
    rf_output_path = os.path.join(rf_folder, name)
    if not isinstance(tree_object_a, Tree):
        tree_object_a = generate_tree_object_from_newick(tree_object_a)
    if not isinstance(tree_object_b, Tree):
        tree_object_b = generate_tree_object_from_newick(tree_object_b)
    rf_first_phase_trees = unify_text_files([tree_object_a.write(format=1), tree_object_b.write(format=1)],
                                            rf_output_path, str_given=True)
    rf = calculate_rf_dist(rf_first_phase_trees, rf_folder,
                           prefix=name)
    return rf


def extract_raxml_statistics_from_msa(full_file_path, output_name, msa_stats, curr_run_directory):
    check_validity_prefix = os.path.join(curr_run_directory, output_name + "_CHECK")
    check_validity_command = (
        "{raxml_exe_path} {threads_config} --force msa --force perf_threads --check --msa {msa_path} --model {model}+G --prefix {prefix} --seed {seed}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_ng_command_prefix(),
        msa_path=full_file_path, prefix=check_validity_prefix, seed=SEED, model = msa_stats["evo_model"])
    reduced_file = check_validity_prefix + ".raxml.reduced.phy"
    # check_log_path = check_validity_prefix + ".raxml.log"
    execute_commnand_and_write_to_log(check_validity_command, run_locally=True)
    if os.path.exists(reduced_file):
        logging.error("Need to re-calculate data on reduced version in " + reduced_file)
        msa_stats["orig_reduced_file_path"] = reduced_file
        raise RE_RUN_ON_REDUCED_VERSION("Input MSA is not valid, re-running on a reduced version")
    parse_prefix = os.path.join(curr_run_directory, output_name + "_PARSE")
    parse_command = "{raxml_exe_path} --force msa --force perf_threads {threads_config} --parse --msa {msa_path} --model {model}+G --prefix {prefix} --seed {seed}".format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_ng_command_prefix(),
        msa_path=full_file_path, prefix=parse_prefix, seed=SEED,model = msa_stats["evo_model"])
    execute_commnand_and_write_to_log(parse_command, run_locally=True)
    binary_msa = parse_prefix + ".raxml.rba"
    msa_stats["local_binary_msa"] = binary_msa
    parsimony_tree_generation_prefix = os.path.join(curr_run_directory, output_name + "pars")
    constant_branch_length_parsimony_tree_path = parsimony_tree_generation_prefix + ".raxml.startTree"
    parsimony_tree_generation_command = (
        "{raxml_exe_path} {threads_config} --force msa --force perf_threads --start --msa {msa_path} --model {model}+G --tree pars{{{n_parsimony_trees}}} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_ng_command_prefix(),
        msa_path=full_file_path, n_parsimony_trees=1, prefix=parsimony_tree_generation_prefix, seed=SEED,model = msa_stats["evo_model"])
    execute_commnand_and_write_to_log(parsimony_tree_generation_command, run_locally=True)
    msa_stats["raxml_parsimony_tree_path"] = constant_branch_length_parsimony_tree_path
    wait_for_file_existence(constant_branch_length_parsimony_tree_path, "Parsimony tree")
    parsimony_model_evaluation_prefix = os.path.join(curr_run_directory, output_name + "pars_eval")
    parsimony_model_and_bl_evaluation_command = (
        "{raxml_exe_path} {threads_config} --force msa --force perf_threads --evaluate --msa {msa_path} --model {model}+G  --tree {parsimony_tree_path} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_ng_command_prefix(),
        msa_path=full_file_path, parsimony_tree_path=constant_branch_length_parsimony_tree_path, seed=SEED,
        prefix=parsimony_model_evaluation_prefix,model = msa_stats["evo_model"])
    execute_commnand_and_write_to_log(parsimony_model_and_bl_evaluation_command, run_locally=True)
    parsimony_log_path = parsimony_model_evaluation_prefix + ".raxml.log"
    wait_for_file_existence(parsimony_log_path, "Parsimony log")
    parsimony_optimized_tree_path = parsimony_model_evaluation_prefix + ".raxml.bestTree"
    parsimony_optimized_model = parsimony_model_evaluation_prefix + ".raxml.bestModel"
    msa_stats["parsimony_optimized_tree_path"] = parsimony_optimized_tree_path
    parsimony_divergence = compute_tree_divergence(parsimony_optimized_tree_path)
    parsimony_tree_alpha = extract_param_from_raxmlNG_log(parsimony_log_path, "alpha")
    if not msa_stats["use_raxml_search"]:
        mad_command = "{mad_exe_path} -t -s {tree_path}".format(mad_exe_path=MAD_COMMAND_PREFIX,
                                                                tree_path=parsimony_optimized_tree_path)
        execute_commnand_and_write_to_log(mad_command, run_locally=True)
        mad_log_path = parsimony_optimized_tree_path + ".rooted"
        mad = extract_mad_file_statistic(mad_log_path)
        msa_stats["mad"] = mad
    msa_stats["alpha"] = parsimony_tree_alpha
    msa_stats["divergence"] = parsimony_divergence
    msa_stats["pars_optimized_model"] = parsimony_optimized_model


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

def generate_n_tree_neighbours_topology_optimized_brlen(n, alpha, original_file_path, curr_run_directory,
                                                   curr_msa_stats, seed):
    prefix = os.path.join(curr_run_directory, "rand_top")
    if curr_msa_stats["use_parsimony_training_trees"]:
        tree_type = "pars"
    else:
        tree_type = "rand"
    random_tree_generation_command = (
        "{raxml_exe_path} {threads_config} --force msa --force perf_threads  --msa {msa_path} --model {model} --start --tree {tree_type}{{{n}}} --prefix {prefix} --opt-branches off --seed {seed} ").format(
        n=1, raxml_exe_path=RAXML_NG_EXE, tree_type=tree_type,
        threads_config=generate_raxml_ng_command_prefix(cpus=curr_msa_stats["n_cpus_training"]),
        msa_path=original_file_path, alpha=alpha, prefix=prefix, seed=seed, model = curr_msa_stats["pars_optimized_model"])
    random_tree_path = prefix + ".raxml.startTree"
    raxml_log_file = prefix + ".raxml.log"
    execute_commnand_and_write_to_log(random_tree_generation_command, curr_run_directory,
                                      job_folder_name="generate_random_trees_job",
                                      job_name="rand_trees", log_file_path=raxml_log_file,
                                      cpus=1, nodes=1, queue=curr_msa_stats["queue"],
                                      run_locally=curr_msa_stats["run_raxml_commands_locally"])
    wait_for_file_existence(random_tree_path, "random tree")
    elapsed_running_time = extract_param_from_raxmlNG_log(raxml_log_file, 'time')
    random_tree_object = generate_tree_object_from_newick(tree_path=random_tree_path)



    if curr_msa_stats["use_parsimony_training_trees"] and n>1:
        logging.debug("Removing duplicates parismony topologies")
        rf_prefix = os.path.join(curr_run_directory, "parsimony_rf_eval")
        rf_command = (
            "{raxml_exe_path} --force msa --force perf_threads --rfdist --tree {rf_file_path} --prefix {prefix}").format(
            raxml_exe_path=RAXML_NG_EXE, rf_file_path=random_tree_path, prefix=rf_prefix)
        execute_commnand_and_write_to_log(rf_command, run_locally=True)
        rf_distances_file_path = rf_prefix + ".raxml.rfDistances"
        random_tree_path = extract_parsimony_unique_topologies(curr_run_directory, random_tree_path,
                                                               rf_distances_file_path, n)
    return random_tree_path, elapsed_running_time

def generate_n_random_tree_topology_constant_brlen(n, alpha, original_file_path, curr_run_directory,
                                                   curr_msa_stats, seed):
    prefix = os.path.join(curr_run_directory, "rand_top")
    if curr_msa_stats["use_parsimony_training_trees"]:
        tree_type = "pars"
    else:
        tree_type = "rand"
    random_tree_generation_command = (
        "{raxml_exe_path} {threads_config} --force msa --force perf_threads  --msa {msa_path} --model {model} --start --tree {tree_type}{{{n}}} --prefix {prefix} --opt-branches off --seed {seed} ").format(
        n=n, raxml_exe_path=RAXML_NG_EXE, tree_type=tree_type,
        threads_config=generate_raxml_ng_command_prefix(cpus=curr_msa_stats["n_cpus_training"]),
        msa_path=original_file_path, alpha=alpha, prefix=prefix, seed=seed, model = curr_msa_stats["pars_optimized_model"])
    random_tree_path = prefix + ".raxml.startTree"
    raxml_log_file = prefix + ".raxml.log"
    execute_commnand_and_write_to_log(random_tree_generation_command, curr_run_directory,
                                      job_folder_name="generate_random_trees_job",
                                      job_name="rand_trees", log_file_path=raxml_log_file,
                                      cpus=1, nodes=1, queue=curr_msa_stats["queue"],
                                      run_locally=curr_msa_stats["run_raxml_commands_locally"])
    wait_for_file_existence(random_tree_path, "random tree")
    elapsed_running_time = extract_param_from_raxmlNG_log(raxml_log_file, 'time')
    if curr_msa_stats["use_parsimony_training_trees"] and n>1:
        logging.debug("Removing duplicates parismony topologies")
        rf_prefix = os.path.join(curr_run_directory, "parsimony_rf_eval")
        rf_command = (
            "{raxml_exe_path} --force msa --force perf_threads --rfdist --tree {rf_file_path} --prefix {prefix}").format(
            raxml_exe_path=RAXML_NG_EXE, rf_file_path=random_tree_path, prefix=rf_prefix)
        execute_commnand_and_write_to_log(rf_command, run_locally=True)
        rf_distances_file_path = rf_prefix + ".raxml.rfDistances"
        random_tree_path = extract_parsimony_unique_topologies(curr_run_directory, random_tree_path,
                                                               rf_distances_file_path, n)
    return random_tree_path, elapsed_running_time


def extract_parsimony_unique_topologies(curr_run_directory, trees_path, dist_path, n):
    rf_prefix = os.path.join(curr_run_directory, "parsimony_rf")
    rf_command = (
        "{raxml_exe_path} --force msa --force perf_threads --rfdist --tree {rf_file_path} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE, rf_file_path=trees_path, prefix=rf_prefix)
    execute_commnand_and_write_to_log(rf_command, run_locally=True)
    unique_file_path = trees_path + "_unique"
    unique_topology_inds = set(list(range(n)))
    with open(dist_path, 'r') as DIST, open(trees_path, 'r') as TREES, open(unique_file_path, 'w') as UNIQUE_TREES:
        distances = DIST.readlines()
        original_trees = TREES.readlines()
        for line in distances:
            lst = line.split("\t")
            curr_tree, comp_tree, dist = int(lst[0]), int(lst[1]), int(lst[2])
            if curr_tree in unique_topology_inds and comp_tree in unique_topology_inds and dist == 0:
                unique_topology_inds.remove(comp_tree)
        unique_trees = [original_trees[ind] for ind in unique_topology_inds]
        n_unique_top = len(unique_trees)
        logging.debug(f'Found {n_unique_top} unique topologies')
        UNIQUE_TREES.writelines(unique_trees)
    rf_prefix = os.path.join(curr_run_directory, "parsimony_check_rf")
    rf_command = (
        "{raxml_exe_path} --force msa --force perf_threads --rfdist --tree {rf_file_path} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE, rf_file_path=unique_file_path, prefix=rf_prefix)
    execute_commnand_and_write_to_log(rf_command, run_locally=True)
    return unique_file_path


def filter_unique_topologies(curr_run_directory, trees_path, n):
    logging.debug("Removing duplicate SPR neighbours")
    rf_prefix = os.path.join(curr_run_directory, "SPR_neighbours")
    rf_command = (
        "{raxml_exe_path} --force msa --force perf_threads --rfdist --tree {rf_file_path} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE, rf_file_path=trees_path, prefix=rf_prefix)
    execute_commnand_and_write_to_log(rf_command, run_locally=True)
    rf_distances_file_path = rf_prefix + ".raxml.rfDistances"
    unique_file_path = trees_path + "_unique"
    unique_topology_inds = set(list(range(n)))
    with open(rf_distances_file_path, 'r') as DIST, open(trees_path, 'r') as TREES, open(unique_file_path,
                                                                                         'w') as UNIQUE_TREES:
        distances = DIST.readlines()
        original_trees = TREES.readlines()
        for line in distances:
            lst = line.split("\t")
            curr_tree, comp_tree, dist = int(lst[0]), int(lst[1]), int(lst[2])
            if curr_tree in unique_topology_inds and comp_tree in unique_topology_inds and dist == 0:
                unique_topology_inds.remove(comp_tree)
        unique_trees = [original_trees[ind] for ind in unique_topology_inds]
        n_unique_top = len(unique_trees)
        logging.debug(f'Found {n_unique_top} unique topologies')
        UNIQUE_TREES.writelines(unique_trees)
    rf_prefix = os.path.join(curr_run_directory, "SPR_neighbours_check")
    rf_command = (
        "{raxml_exe_path} --force msa --force perf_threads --rfdist --tree {rf_file_path} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE, rf_file_path=unique_file_path, prefix=rf_prefix)
    execute_commnand_and_write_to_log(rf_command, run_locally=True)
    return unique_file_path


def raxml_compute_tree_per_site_ll(curr_run_directory, full_data_path, tree_file, ll_on_data_prefix, alpha,
                                   curr_msa_stats,
                                   opt_brlen=True, model = None):
    compute_site_ll_prefix = os.path.join(curr_run_directory, ll_on_data_prefix)
    if not model:
        model = curr_msa_stats["pars_optimized_model"]
    brlen_command = "--opt-branches off --opt-model off " if not opt_brlen else ""
    compute_site_ll_run_command = (
        "{raxml_exe_path} {threads_config} --force msa --force perf_threads --sitelh --msa {msa_path} --model {model} {brlen_command} --tree {tree_file} --seed {seed} --prefix {compute_site_ll_prefix} ").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_ng_command_prefix(cpus=curr_msa_stats["n_cpus_training"]),
        alpha=alpha, msa_path=full_data_path, tree_file=tree_file, seed=SEED,
        prefix=compute_site_ll_prefix, brlen_command=brlen_command, compute_site_ll_prefix=compute_site_ll_prefix, model = curr_msa_stats["pars_optimized_model"])
    sitelh_file = compute_site_ll_prefix + ".raxml.siteLH"
    raxml_log_file = compute_site_ll_prefix + ".raxml.log"
    execute_commnand_and_write_to_log(compute_site_ll_run_command, curr_run_directory,
                                      job_folder_name="raxml_ll_eval_job_for_training",
                                      job_name="training_opt", log_file_path=raxml_log_file,
                                      cpus=curr_msa_stats["n_cpus_training"], nodes=curr_msa_stats["n_nodes_training"],
                                      queue=curr_msa_stats["queue"],
                                      run_locally=curr_msa_stats["run_raxml_commands_locally"])
    wait_for_file_existence(sitelh_file, "Sitelh file")
    sitelh_list = raxml_extract_sitelh(sitelh_file)
    elapsed_running_time = extract_param_from_raxmlNG_log(raxml_log_file, 'time')
    return sitelh_list, elapsed_running_time


def raxml_optimize_trees_for_given_msa(full_data_path, ll_on_data_prefix, tree_file, msa_stats,
                                       curr_run_directory, opt_brlen=True, weights=None, return_trees_file=False,
                                       n_cpus=1, model = None):
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
    brlen_command = "--opt-branches off --opt-model off " if not opt_brlen else ""
    if model is None:
        model = f"{msa_stats['evo_model']}+G{{{msa_stats['alpha']}}}" if not msa_stats.get("pars_optimized_model") else msa_stats["pars_optimized_model"]
    compute_ll_run_command = (
        "{raxml_exe_path} {threads_config} --force msa --force perf_threads --evaluate --msa {msa_path} --model {model} {brlen_command} --tree {tree_file} {weights_path_command} --seed {seed} --prefix {prefix}").format(
        raxml_exe_path=RAXML_NG_EXE,
        threads_config=generate_raxml_ng_command_prefix(n_cpus),
        alpha=alpha, msa_path=full_data_path, tree_file=tree_file, seed=SEED,
        prefix=prefix, weights_path_command=weights_path_command, brlen_command=brlen_command, model = model)
    optimized_trees_path = prefix + ".raxml.mlTrees"
    best_tree_path = prefix + ".raxml.bestTree"
    raxml_log_file = prefix + ".raxml.log"
    execute_commnand_and_write_to_log(compute_ll_run_command, curr_run_directory,
                                      job_folder_name="raxml_optimize_test_trees_job",
                                      job_name="trees_opt", log_file_path=raxml_log_file,
                                      cpus=n_cpus, nodes=msa_stats["n_nodes_training"],
                                      queue=msa_stats["queue"], run_locally=msa_stats["run_raxml_commands_locally"])

    trees_ll_on_data = extract_param_from_raxmlNG_log(raxml_log_file, "ll")
    elapsed_running_time = extract_param_from_raxmlNG_log(raxml_log_file, 'time')
    optimized_trees_final_path = optimized_trees_path if os.path.exists(optimized_trees_path) else best_tree_path
    tree_objects = generate_multiple_tree_object_from_newick(optimized_trees_final_path)
    if return_trees_file:
        return trees_ll_on_data, optimized_trees_final_path, elapsed_running_time
    return trees_ll_on_data, tree_objects, elapsed_running_time





