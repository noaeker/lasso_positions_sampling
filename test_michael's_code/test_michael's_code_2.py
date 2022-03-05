from os import path
import sys
import os

from ete3 import *
from spr_prune_and_regraft import *
from raxml import *
from help_functions import *
from generate_SPR import *
import re

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



def str_to_spr_move(edge_str):
    match = re.search(r'\[a=(\w+) b=(\w+)\],\[a=(\w+) b=(\w+)\]',edge_str)
    a1,b1,a2,b2 = match.group(1),match.group(2),match.group(3),match.group(4)
    return (Edge(node_a=a1,node_b=b1), Edge(node_a=a2,node_b=b2))

def main():
    tree_orig = Tree('(Sp002:0.00121364,Sp008:0.00376252,((Sp005:1e-08,Sp004:1e-08)N4:0.00246211,(Sp009:0.00366749,((Sp001:0.0432468,(Sp006:1e-08,(Sp007:0.00117634,Sp003:2.5e-09)N15:2.5e-09)N13:0.00245986)N10:0.00245986,Sp000:0.00862005)N9:0.00626271)N5:0.00626271)N3:0.0457089);', format=1)
    tree_orig.get_tree_root().name = "ROOT"
    msa_path = "masked_species_real_msa.phy"
    print( tree_orig.get_ascii(attributes=['name'], show_internal=True) )
    orig_tree_path = write_tree_objects_to_file(tree_orig,
                                                os.getcwd(),
                                                "tree_orig")
    orig_data = pd.read_csv("two_step_lookahead_QA.csv")
    first_moves = list(orig_data["first_move"])
    second_moves = list(orig_data["second_move"])
    for first_move, second_move in zip(first_moves,second_moves):
        first_move = str_to_spr_move(first_move)
        neighbor = generate_neighbour(tree_orig, first_move)
        print(neighbor.get_ascii(attributes=['name'], show_internal=True))
        print(second_move)
        second_move = str_to_spr_move(second_move)
        neighbor2 = generate_neighbour(neighbor, second_move)
        #neighbour = generate_neighbour(tree_orig, possible_move = (Edge(node_a="N1", node_b="ROOT"), Edge(node_a="N8", node_b="N6")) )






if __name__ == "__main__":
    main()





