from os import path
import sys
import os

from ete3 import *
from spr_prune_and_regraft import *
from raxml import *
from help_functions import *
from generate_SPR import *

RAXML_NG_EXE = "/Users/noa/Programs/Raxml/raxml-ng  "

def RAxML_NG_EVAL(msa_results_directory, raxml_ng_exe_path, msa_path, tree_path,name):
    '''
    RAxML-NG is run with default parameters on msa_path, and final tree log-likelihood is evaluated with respect to true_msa_path
    '''
    prefix = os.path.join(msa_results_directory, name)
    ll_prefix = os.path.join(msa_results_directory, f"true_ll_{name}")
    compute_ll_run_command = (
    f"{raxml_ng_exe_path} --evaluate --msa {msa_path} --model GTR+G --tree {tree_path} --seed 1 --prefix {ll_prefix} ")
    subprocess.run(compute_ll_run_command, shell=True)
    true_ll_log_file = ll_prefix + ".raxml.log"
    with open(true_ll_log_file) as TRUE_LL_LOG_FILE:
        data = TRUE_LL_LOG_FILE.read()
    pattern = r'Tree #\d+, final logLikelihood: (-[\d.]+)'
    ll_string = re.findall(pattern, data)[0]
    ll_float = float(ll_string)
    return ll_float


def main():
    tree_orig = Tree('((((((Sp005:0.00322175,((Sp011:1e-08,Sp000:0.00238396)N18:0.00079083,(Sp010:0.00129242,Sp003:0.00255366)N19:0.00063641)N17:5e-08)N12:4e-08,Sp002:0.0170639)N8:0.00195312,Sp009:0.0206418)N6:0.00354918,((Sp012:1e-08,Sp004:4e-08)N10:4e-08,Sp007:0.0133045)N7:0.0136031)N4:1.5e-07,Sp008:0.0524783)N1:0.0827905,Sp006:0.00716386,Sp001:0.117007);', format=1)
    tree_orig.get_tree_root().name = "ROOT"
    msa_path = "masked_species_real_msa.phy"
    print( tree_orig.get_ascii(attributes=['name'], show_internal=True) )
    orig_tree_path = write_tree_objects_to_file(tree_orig,
                                                os.getcwd(),
                                                "tree_orig")

    tree_greedy1 = generate_neighbour(tree_orig, possible_move = (Edge(node_a="N1", node_b="ROOT"), Edge(node_a="N8", node_b="N6")) )
    #print(tree_greedy1.get_ascii(attributes=['name'], show_internal=True))
    greedy1_tree_path = write_tree_objects_to_file(tree_greedy1,
                                                os.getcwd(),
                                                "tree_greedy1")

    tree_greedy2 = generate_neighbour(tree_greedy1,
                                      possible_move=(Edge(node_a="Sp002", node_b="N8"), Edge(node_a="Sp003", node_b="N19")))
    #print(tree_greedy2.get_ascii(attributes=['name'], show_internal=True))
    greedy2_tree_path = write_tree_objects_to_file(tree_greedy2,
                                                   os.getcwd(),
                                                   "tree_greedy2")

    tree_not_greedy1 = generate_neighbour(tree_orig,
                                      possible_move=(Edge(node_a="N12", node_b="N8"), Edge(node_a="N1", node_b="ROOT")))
    #print( tree_not_greedy1 .get_ascii(attributes=['name'], show_internal=True))
    tree_not_greedy1_path = write_tree_objects_to_file(tree_not_greedy1,
                                                   os.getcwd(),
                                                   "tree_not_greedy1")

    tree_not_greedy2 = generate_neighbour(tree_not_greedy1,
                                      possible_move=(Edge(node_a="N6", node_b="N4"), Edge(node_a="N18", node_b="N17")))
    #print( tree_not_greedy2 .get_ascii(attributes=['name'], show_internal=True))
    tree_not_greedy2_path = write_tree_objects_to_file(tree_not_greedy2,
                                                   os.getcwd(),
                                                   "tree_not_greedy2")
    #tree_greedy2 = generate_neighbour(tree_orig, possible_move = (Edge(node_a="N3", node_b="ROOT"), Edge(node_a="N4", node_b="N3")) )
    #tree_non_greedy1 = generate_neighbour(tree_orig, possible_move = (Edge(node_a="N3", node_b="ROOT"), Edge(node_a="N4", node_b="N3")) )
    #tree_non_greedy2 = generate_neighbour(tree_orig, possible_move = (Edge(node_a="N3", node_b="ROOT"), Edge(node_a="N4", node_b="N3")) )


    results_folder = os.path.join(os.getcwd(),"raxml_outputs")
    create_or_clean_dir(results_folder)
    orig_tree_ll = RAxML_NG_EVAL(results_folder, RAXML_NG_EXE, msa_path, tree_path = orig_tree_path,name = "orig_tree_eval")
    greedy1_tree_ll = RAxML_NG_EVAL(results_folder, RAXML_NG_EXE, msa_path, tree_path=greedy1_tree_path,
                                 name="greedy1_tree_eval")
    greedy2_tree_ll = RAxML_NG_EVAL(results_folder, RAXML_NG_EXE, msa_path, tree_path=greedy2_tree_path,
                                    name="greedy2_tree_eval")
    not_greedy1_tree_ll = RAxML_NG_EVAL(results_folder, RAXML_NG_EXE, msa_path, tree_path=tree_not_greedy1_path,
                                    name="not_greedy1_tree_eval")
    not_greedy2_tree_ll = RAxML_NG_EVAL(results_folder, RAXML_NG_EXE, msa_path, tree_path=tree_not_greedy2_path,
                                    name="not_greedy2_tree_eval")

    print(f"orig_tree_ll = {orig_tree_ll}")
    print(f"greedy1_tree_ll = {greedy1_tree_ll}")
    print(f"greedy2_tree_ll = {greedy2_tree_ll}")
    print(f"Not greedy1_tree_ll = {not_greedy1_tree_ll}")
    print(f"Not greedy2_tree_ll = {not_greedy2_tree_ll}")




if __name__ == "__main__":
    main()





