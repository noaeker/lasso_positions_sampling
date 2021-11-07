
from help_functions import *
import re
import os

def get_rate4site(MSA_path, tree_path, output_path):
    rate4site_command = "{rate4site_exe_path}  -bn -s {msa_path} -t {tree_path} -o {output_path}".format(
        rate4site_exe_path = RATE4SITE_COMMAND_PREFIX, msa_path = MSA_path, tree_path = tree_path, output_path = output_path)
    subprocess.run(rate4site_command, shell=True)
    if os.path.exists("r4s.res"):
        os.remove("r4s.res")
    if os.path.exists("r4sOrig.res"):
        os.remove("r4sOrig.res")
    return rate4site_command


def parse_rate4site(rate4site_output_path):
    scores = []
    with open(rate4site_output_path,'r') as RATE4SITE_OUTPUT_PATH:
        rate4site = RATE4SITE_OUTPUT_PATH.readlines()
        for line in rate4site:
            match = re.search('\d+\s+.\s+([\d.e-]+\s+)',line)
            if match:
                scores.append(float(match.group(1)))
    return scores




