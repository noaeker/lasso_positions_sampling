import pickle
from raxml import *
from help_functions import *
import pandas as pd

NJOBS=1
#OPT_BRLEN = True
#brlen_command = "--opt-branches off" if not OPT_BRLEN else ""
PREFIX  = "/Users/noa/Workspace/lasso_positions_sampling_results/test_diff_change_2"
CHOSEN_TRAINING_SIZE=800
LL_CALCULATION_COMMAND = "{raxml_exe_path} --evaluate   --msa {msa_path} --model WAG+G{{{alpha}}} --tree {tree_file} {brlen_command} {weights_path_command} --prefix {prefix}"
CONST = 10000000
CHOSEN_BRLEN = "optimized"


def get_ll_and_tree(prefix,tree_file,msa_path,brlen_command="", weight_path_command=""):
    command = LL_CALCULATION_COMMAND.format(
        raxml_exe_path=RAXML_NG_EXE,
        alpha=alpha, msa_path=msa_path, tree_file=tree_file,
        prefix=prefix, brlen_command=brlen_command,weights_path_command=weight_path_command)
    subprocess.run(command, shell=True)
    raxml_log_file = prefix + ".raxml.log"
    tree_path = prefix + ".raxml.bestTree"
    ll = extract_param_from_raxmlNG_log(raxml_log_file, "ll")
    return ll,tree_path


for job_ind in range(NJOBS):
    job_dict={}

    job_dir=(os.path.join(PREFIX,"job_{}".format(job_ind)))
    internal_dir=([name for name in os.listdir(job_dir) if "Selectome" in name][0])
    running_dir = os.path.join(job_dir,internal_dir)
    sampled_file_path = '{job_folder}/Lasso_folder/{brlen}/training_{training_size}_random_tree_eval/0_fixed_sampled.phy'.format(job_folder=running_dir, training_size= CHOSEN_TRAINING_SIZE, brlen=CHOSEN_BRLEN)
    lasso_weights_file_path = '{job_folder}/Lasso_folder/{brlen}/training_{training_size}_random_tree_eval/0_fixed_lasso_weights.txt'.format(job_folder=running_dir,training_size= CHOSEN_TRAINING_SIZE, brlen= CHOSEN_BRLEN)
    full_msa_path = '{job_folder}/raxml_full_data_results_0/0_fixed.phy'.format(job_folder= running_dir)
    test_tree_path = '{job_folder}/raxml_full_data_results_0_fixed/0pars_eval.raxml.bestTree'.format(job_folder=running_dir)
    curr_results_folder = '{job_dir}/test_issues'.format(job_dir=job_dir)
    create_or_clean_dir(curr_results_folder)
    lasso_model_dump_path = '{job_folder}/Lasso_folder/{brlen}/training_{training_size}_random_tree_eval/lasso_model.sav'.format(job_folder=running_dir,training_size= CHOSEN_TRAINING_SIZE, brlen=CHOSEN_BRLEN)
    weights_file_path_alternative = os.path.join(curr_results_folder,"weights.alternative")
    constant_weight_file_path =  os.path.join(curr_results_folder,"weights.constant")
    msa_stats_dump_path = '{job_folder}/curr_msa_stats.dump'.format(job_folder=running_dir)
    lasso_weights_path_command = "--site-weights {}".format(lasso_weights_file_path)
    constant_weights_path_command = "--site-weights {}".format(constant_weight_file_path)
    weights_path_command_alternative = "--site-weights {}".format(weights_file_path_alternative)

    with open(lasso_model_dump_path, 'rb') as model:
        lasso_model = pickle.load(model)
    intercept=lasso_model.intercept_
    if USE_INTEGER_WEIGHTS:
        weights = [int(lasso_model.coef_[ind] * INTEGER_CONST) for ind in range(len(lasso_model.coef_))]
    else:
        weights = lasso_model.coef_
    lasso_chosen_locis = [ind for ind in range(len(weights)) if weights[ind] != 0]
    with open(lasso_weights_file_path,'r') as weights_f:
        f= weights_f.read().split(" ")
        lasso_weights= [int(n) for n in f if len(n)>0]
    msa_length = len(lasso_model.coef_)
    with open(msa_stats_dump_path,'rb') as model:
        curr_msa_stats = pickle.load(model)
    with open(weights_file_path_alternative, 'w') as f:
       for i in range(msa_length):
           if i in lasso_chosen_locis:
               f.write(str(CONST) + " ")
           else:
               f.write(str(1)+" ")
    with open(constant_weight_file_path, 'w') as f:
        for i in range(msa_length):
            f.write(str(CONST) + " ")

    alpha = curr_msa_stats["alpha"]


    take_half_msa_path = '{curr_results_folder}/half_sampled.phy'.format(curr_results_folder=curr_results_folder)
    write_to_sampled_alignment_path(curr_msa_stats["alignment_data"], take_half_msa_path, range(curr_msa_stats["n_loci"] // 2),
                                    curr_msa_stats["file_type_biopython"])
    take_half_weights_path = '{curr_results_folder}/half_sampled_weights.phy'.format(curr_results_folder=curr_results_folder)
    with open(take_half_weights_path, 'w') as f:
        for i in range(msa_length//2):
            f.write(str(CONST) + " ")
        for i in range(msa_length//2,msa_length):
            f.write("1"+" ")

    half_weights_path_command = "--site-weights {}".format(take_half_weights_path)

    prefix_standard_opt =  os.path.join(curr_results_folder, "standard_full_data_opt")
    prefix_standard_eval = os.path.join(curr_results_folder, "standard_full_data_eval")
    prefix_constant_weights = os.path.join(curr_results_folder, "constant_weights_full_data")
    prefix_alternative = os.path.join(curr_results_folder, "alternative_weights_full_data")
    prefix_alternative_eval = os.path.join(curr_results_folder, "alternative_weights_eval")
    prefix_lasso =  os.path.join(curr_results_folder, "lasso_weights_sampled_data")
    prefix_standard_on_lasso = os.path.join(curr_results_folder, "eval_lasso_tree_using_standard")
    prefix_lasso_on_standard = os.path.join(curr_results_folder, "eval_standard_tree_using_lasso")
    prefix_lasso_on_lasso = os.path.join(curr_results_folder, "eval_lasso_tree_using_lasso")
    prefix_sampled_data = os.path.join(curr_results_folder, "opt_sampled_data")
    prefix_half_sampled_data = os.path.join(curr_results_folder, "opt_half_sampled_data")
    prefix_half_weights_data = os.path.join(curr_results_folder, "opt_half_weights_data")

    job_dict["job_ind"] = job_ind
    job_dict["test_tree_raw"] = get_tree_string(test_tree_path)
    job_dict["standard_ll_OPT"],standard_tree_OPT=get_ll_and_tree(prefix_standard_opt, test_tree_path, full_msa_path)
    job_dict["standard_tree_OPT_newick"] = get_tree_string(standard_tree_OPT)
    job_dict["standard_ll_EVAL"], standard_tree_EVAL_path = get_ll_and_tree(prefix_standard_eval,
                                                                            test_tree_path, full_msa_path, brlen_command= "--opt-branches off")

    job_dict["standard_tree_EVAL_newick"] = get_tree_string(standard_tree_EVAL_path)
    job_dict["optimized_standard_tree_EVAL_by_LASSO"] = get_ll_and_tree(prefix_lasso_on_standard,
                                                                                   standard_tree_OPT, sampled_file_path,
                                                                                   brlen_command="--opt-branches off", weight_path_command=lasso_weights_path_command)[0]/INTEGER_CONST+intercept

    alternative_weights_ll_OPT, alternative_weights_tree_OPT = get_ll_and_tree(prefix_alternative,
                                                                        test_tree_path,full_msa_path,
                                                                        brlen_command="",
                                                                        weight_path_command=weights_path_command_alternative)
    job_dict["alternative_weights_ll_OPT"] = alternative_weights_ll_OPT / CONST
    job_dict["alternative_weights_OPT_newick"] = get_tree_string(alternative_weights_tree_OPT)

    constant_weights_ll_OPT, constant_weights_tree_OPT = get_ll_and_tree(prefix_constant_weights,
                                                                                                       test_tree_path,full_msa_path,
                                                                                                       brlen_command="",
                                                                                                       weight_path_command=constant_weights_path_command)

    job_dict["constant_weights_ll_OPT"] = constant_weights_ll_OPT / CONST
    job_dict["constant_weights_OPT_newick"] = get_tree_string(constant_weights_tree_OPT)
    job_dict["sampled_data_ll_OPT"], sampled_data_tree_OPT = get_ll_and_tree(prefix_sampled_data ,test_tree_path,
                                                                                                       sampled_file_path
                                                                                                    )

    job_dict["sampled_data_OPT_newick"] = get_tree_string(sampled_data_tree_OPT)

    # job_dict["half_sampled_data_ll_OPT"], half_sampled_data_tree_OPT = get_ll_and_tree(prefix_half_sampled_data ,test_tree_path,
    #                                                                                                   take_half_msa_path
    #                                                                                                 )
    #
    # job_dict["half_sampled_data_OPT_newick"] = get_tree_string(half_sampled_data_tree_OPT )
    #
    # job_dict["half_weights_opt_ll"], half_weight_tree_OPT = get_ll_and_tree(  prefix_half_weights_data ,
    #                                                                                    test_tree_path,
    #                                                                                    full_msa_path,
    #                                                                           weight_path_command= half_weights_path_command
    #                                                                                    )
    #
    # job_dict["half_weight_data_OPT_newick"] = get_tree_string(half_weight_tree_OPT )


    sitelh_list = raxml_compute_tree_per_site_ll(curr_results_folder, full_msa_path, test_tree_path,
                                                 "standard_OPT_sitelh", alpha, opt_brlen=True)

    job_dict["ll_Lasso_using_standard_sitelh"] = (sum([lasso_weights[i]*sitelh_list[j] for i,j in enumerate(lasso_chosen_locis)])/INTEGER_CONST)+intercept


    Lasso_ll_OPT, Lasso_tree_OPT = get_ll_and_tree(prefix_lasso,test_tree_path,
                                                                                                       sampled_file_path,brlen_command="",
                                                                           weight_path_command=lasso_weights_path_command)
    job_dict["Lasso_ll_OPT"] = Lasso_ll_OPT/INTEGER_CONST+intercept


    job_dict["Lasso_OPT_newick"] = get_tree_string( Lasso_tree_OPT)

    job_dict["Lasso_tree_EVAL_by_standard"] = get_ll_and_tree(prefix_standard_on_lasso, Lasso_tree_OPT,
                                                                           full_msa_path, brlen_command="--opt-branches off"
                                                                           )[0]
    job_dict["Lasso_tree_EVAL_by_Lasso"] = get_ll_and_tree(prefix_lasso_on_lasso, Lasso_tree_OPT,
                                                              sampled_file_path, brlen_command="--opt-branches off",weight_path_command=lasso_weights_path_command
                                                              )[0]/INTEGER_CONST+intercept

    print(job_dict)
