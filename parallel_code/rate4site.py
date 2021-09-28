
from help_functions import *
import re
import os

def get_rate4site(MSA_path, tree_path, output_path):
    rate4site_command = "{rate4site_exe_path}  -bn -s {msa_path} -t {tree_path} -o {output_path}".format(
        rate4site_exe_path = RATE4SITE_COMMAND_PREFIX, msa_path = MSA_path, tree_path = tree_path, output_path = output_path)
    subprocess.run(rate4site_command, shell=True)
    os.remove("r4s.res")
    os.remove("r4sOrig.res")


def parse_rate4site(rate4site_output_path):
    with open(rate4site_output_path,'r') as RATE4SITE_OUTPUT_PATH:
        rate4site = RATE4SITE_OUTPUT_PATH.read()
        score_per_position = re.findall('\d+\s+.\s+([\d.-]+\s+)',rate4site)
        scores = [float(ll) for ll in score_per_position]
    return scores




# MSA_path = "/Users/noa/Workspace/lasso_positions_sampling_results/test_curr/job_0/_Users_noa_Workspace_data_supermatrices_edited__supermatrices_NagyA1_NagyA1.fasta/0.fasta"
# tree_path = "/Users/noa/Workspace/lasso_positions_sampling_results/test_curr/job_0/_Users_noa_Workspace_data_supermatrices_edited__supermatrices_NagyA1_NagyA1.fasta/raxml_full_data_results_0/0pars_eval.raxml.bestTree"
# output_path = "/Users/noa/Workspace/lasso_positions_sampling_results/test_curr/job_0/_Users_noa_Workspace_data_supermatrices_edited__supermatrices_NagyA1_NagyA1.fasta/r4s.res"
# get_rate4site(MSA_path, tree_path, output_path)
# scores = parse_rate4site(output_path)
# print(scores)
# print(len(scores))
# scores = parse_rate4site("r4sOrig.res")
# print(scores)
# print(len(scores))