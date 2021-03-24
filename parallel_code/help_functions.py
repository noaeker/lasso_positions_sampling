import os
import pandas as pd
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from config import *
import argparse

def remove_MSAs_with_not_enough_seq(file_path_list,min_seq):
    proper_file_path_list = []
    for path in file_path_list:
        file_type_biopython = extract_file_type(path, True)
        with open(path) as file:
            n_seq = len(list(SeqIO.parse(file, file_type_biopython)))
            if n_seq>=min_seq:
                proper_file_path_list.append(path)
    return proper_file_path_list

def write_to_sampled_alignment_path(original_alignment_data, sampled_alignment_path, samp_indexes, file_type):
    sampled_sequence = []
    for original_record in original_alignment_data:
        sampled_seq = Seq(''.join([original_record.seq[ind] for ind in samp_indexes]))
        sampled_record = SeqRecord(sampled_seq, id=original_record.id, name=original_record.name,
                                   description=original_record.description)
        sampled_sequence.append(sampled_record)
    val = SeqIO.write(sampled_sequence, sampled_alignment_path, file_type)
    if val == len(original_alignment_data):
        logging.info("   #Sampled columns written succesfully to new file " + sampled_alignment_path)
    else:
        logging.error("   #ERROR: Sampled columns not written succesfully to file " + sampled_alignment_path)


def take_up_to_x_sequences(original_alignment_data,trimmed_alignment_path, number_of_sequences,file_type):
    sampled_sequence = []
    seq_values = set()
    random.seed(SEED)
    random.shuffle(original_alignment_data)
    for record in original_alignment_data:
        if len(sampled_sequence)>=number_of_sequences:
            break
        if record.seq in seq_values:
            continue
        else:
            seq_values.add(record.seq)
            sampled_sequence.append(record)
    try:
        val = SeqIO.write(sampled_sequence, trimmed_alignment_path, file_type)
        logging.info(" {} sequences written succesfully to new file {}".format(number_of_sequences,trimmed_alignment_path))
    except:
        logging.error("ERROR! {} sequences NOT written succesfully to new file {}".format(number_of_sequences,trimmed_alignment_path))


def extract_file_type(path, change_format=False, ete=False):
    filename, file_extension = os.path.splitext(path)
    if change_format:
        if file_extension == '.phy':
            file_extension = 'iphylip' if ete == True else 'phylip-relaxed'
        elif file_extension == ".fasta":
            file_extension = 'fasta'
        elif file_extension == ".nex":
            file_extension = 'nexus'
    return file_extension


def delete_file_content(file_path):
    with open(file_path, 'w'):
        pass




def extract_alignment_files_from_dir(dir):
    files_list = []
    if os.path.exists(dir):
        for file in os.listdir(dir):
            if file.endswith(".phy") or file.endswith(".fasta") or file.endswith(".nex"):
                files_list.append(os.path.join(dir, file))
    return files_list


def extract_dir_list_from_csv(dir_list_csv_path):
    df = pd.read_csv(dir_list_csv_path)
    df.sort_values(by='nchars',ascending=False,inplace=True)
    dir_list = [os.path.join(MSAs_FOLDER,path) for path in list(df["path"])]
    logging.debug("Number of paths in original csv = {n_paths}".format(n_paths = len(df.index)))
    return dir_list

def extract_alignment_files_from_general_csv(dir_list_csv_path):
    files_list = []
    logging.debug("Extracting alignments from {}".format(dir_list_csv_path))
    dir_list = extract_dir_list_from_csv(dir_list_csv_path)
    for dir in dir_list:
        if os.path.exists(dir):
            for file in os.listdir(dir):
                if (file.endswith(".phy") or file.endswith(".fasta")):
                    files_list.append(os.path.join(dir, file))
                    break
        else:
            logging.error("Following MSA dir does not exist {dir}".format(dir=dir))
    logging.debug("Overalls number of MSAs found in the given directories is: {nMSAs}".format(nMSAs=len(files_list)))
    return files_list

def alignment_list_to_df(alignment_data):
    alignment_list = [list(alignment_data[i].seq) for i in range(len(alignment_data))]
    loci_num = len(alignment_data[0].seq)
    columns = list(range(0, loci_num))
    original_alignment_df = pd.DataFrame(alignment_list, columns=columns)
    return original_alignment_df


def delete_dir_content(dir_path):
    for filename in os.listdir(dir_path):
        file_path = os.path.join(dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)

        except Exception as e:
            logging.error('Failed to delete %s. Reason: %s' % (file_path, e))
            return False
    return True

def create_or_clean_dir(dir):
    if os.path.exists(dir):
        delete_dir_content(dir)
    else:
        os.mkdir(dir)


def create_dir_if_not_exists(dir):
    if not os.path.exists(dir):
            os.mkdir(dir)


def unify_text_files(input_file_path_list, output_file_path):
    with open(output_file_path, 'w') as outfile:
        for fname in input_file_path_list:
            with open(fname) as infile:
                outfile.write(infile.read())



def add_csvs_content(csvs_path_list, unified_csv_path):
    existing_df = [pd.read_csv(unified_csv_path)] if os.path.exists(unified_csv_path) else []
    existing_df_size = pd.read_csv(unified_csv_path).size if os.path.exists(unified_csv_path) else 0
    logging.info('Existing df size is: {}'.format(existing_df_size))
    non_empty_df = [pd.read_csv(f) for f in csvs_path_list if not pd.read_csv(f).empty]
    combined_df = pd.concat(non_empty_df+ existing_df,sort=False)
    combined_df_size = combined_df.size
    logging.info('Combined df size is: {}'.format(combined_df_size))
    combined_df.to_csv(unified_csv_path,index=False)
    return combined_df

def remove_empty_columns(csv_path):
    if os.path.exists((csv_path)):
        df = pd.read_csv(csv_path)
        df = df.dropna(how='all', axis=1)
        df.to_csv(csv_path, index = False)




def main_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_prefix', action='store', type=str, default=CURR_RUN_PREFIX)
    parser.add_argument('--jobs_prefix', action='store', type=str, default=CURR_JOBS_PREFIX)
    parser.add_argument('--n_MSAs', action='store', type=int, default=N_MSAS)
    parser.add_argument('--n_jobs', action='store', type=int, default=N_JOBS)
    parser.add_argument('--first_msa_ind', action='store', type=int, default=0)
    parser.add_argument('--n_random_starting_trees', action='store', type=int, default=N_RANDOM_STARTING_TREES)
    parser.add_argument('--random_trees_training_size', action='store', type=int, default=RANDOM_TREES_TRAINING_SIZE)
    parser.add_argument('--exp_brlen',action='store_true',default= True)
    parser.add_argument('--uni_brlen', action='store_true', default= False)
    parser.add_argument('--opt_brlen', action='store_true', default= True)
    parser.add_argument('--const_brlen', action='store_true', default=False)
    parser.add_argument('--random_trees_test_size', action='store', type=int, default=RANDOM_TREES_TEST_SIZE)
    parser.add_argument('--max_n_seq', action='store', type=int, default=MAX_N_SEQ)
    parser.add_argument('--min_n_seq', action='store', type=int, default=MIN_N_SEQ)
    parser.add_argument('--only_evaluate_lasso', action='store_true',default=False)
    parser.add_argument('--lasso_baseline_run_prefix',action='store', type=str, default=LASSO_BASELINE)
    parser.add_argument('--spr_baseline_run_prefix', action='store', type=str, default=SPR_BASELINE)
    args = parser.parse_args()
    return args

#job_ind, curr_job_folder, job_msa_paths_file, job_csv_path, spr_log_path, general_log_path, curr_job_status_file, max_n_sequences, n_random_starting_trees,only_evaluate_lasso
def job_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job_ind', action='store', type=int)
    parser.add_argument('--curr_job_folder', action='store', type=str)
    parser.add_argument('--max_n_sequences', action='store', type=int)
    parser.add_argument('--n_random_starting_trees', action='store', type=int)
    parser.add_argument('--random_trees_training_size', action='store', type=int)
    parser.add_argument('--exp_brlen',action='store_true',default=False)
    parser.add_argument('--uni_brlen', action='store_true', default=False)
    parser.add_argument('--opt_brlen', action='store_true', default=False)
    parser.add_argument('--const_brlen', action='store_true', default=False)
    parser.add_argument('--random_trees_test_size', action='store', type=int)
    parser.add_argument('--only_evaluate_lasso', action='store_true',default = False)
    parser.add_argument('--run_prefix', action='store', type=str)
    parser.add_argument('--lasso_baseline_run_prefix', action='store', type=str)
    parser.add_argument('--spr_baseline_run_prefix', action='store', type=str)
    args = parser.parse_args()
    return args

def get_job_related_files_paths(curr_job_folder, job_ind):
    job_status_file = os.path.join(curr_job_folder, str(job_ind) + "_status")
    job_csv_path=os.path.join(curr_job_folder, str(job_ind) + ".csv")
    job_msa_paths_file = os.path.join(curr_job_folder, "file_paths_" + str(job_ind))
    spr_log_path = os.path.join(curr_job_folder, "job_" + str(job_ind) + "_spr_log.log")
    general_log_path = os.path.join(curr_job_folder, "job_" + str(job_ind) + "_general_log.log")
    return {"job_status_file" : job_status_file, "job_csv_path":job_csv_path, "job_msa_paths_file" : job_msa_paths_file, "spr_log_path": spr_log_path , "general_log_path": general_log_path}
