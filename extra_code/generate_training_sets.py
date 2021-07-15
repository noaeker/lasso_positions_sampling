from help_functions import *

PATH = '/Users/noa/Workspace/data/LARGE_FILES'
def fileCount(path, extension):
    count = 0
    for root, dirs, files in os.walk(path):
        count += sum(f.endswith(extension) for f in files)
    return count

def generate_training_sets(path):
    msa_files = []
    for dirpath, subdirs, files in os.walk(path):
        for f in files:
            if f.endswith('phy') or f.endswith('fasta'):
                msa_files.append(os.path.join(dirpath, f))
    return msa_files



training_set_paths = generate_training_sets(PATH)
dst = '/Users/noa/Workspace/data/LARGE_FILES_edited'
os.mkdir(dst)

for path in training_set_paths:
    file_type_biopython = extract_file_type(path, True)
    with open(path) as original:
        reduced_local_alignment_data = list(SeqIO.parse(original, file_type_biopython))
        n_seq = len(reduced_local_alignment_data)
        n_loci = len(reduced_local_alignment_data[0].seq)
    print(n_loci, n_seq)
    if n_loci >=10000 and n_seq>=15:
        new_file_name =path.replace('/Users/noa/Workspace/data/LARGE_FILES',"").replace('/',"_")
        shutil.copy2(path, os.path.join(dst,new_file_name))
