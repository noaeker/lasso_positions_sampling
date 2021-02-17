import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--dir_path', action='store', type=str)
parser.add_argument('--log_path', action='store', type=str)
args = parser.parse_args()
path = args.run_prefix
shutil.rmtree(path)