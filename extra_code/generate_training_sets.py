from help_functions import *

PATH = '/Users/noa/Workspace/data/supermatrices'
def fileCount(path, extension):
    count = 0
    for root, dirs, files in os.walk(path):
        count += sum(f.endswith(extension) for f in files)
    return count

def generate_training_sets(path):
    msa_files = []
    for dirpath, subdirs, files in os.walk(path):
        for f in files:
            if f.endswith('phy') or f.endswith('fasta') or f.endswith('aln'):
                msa_files.append(os.path.join(dirpath, f))
    return msa_files







training_set_paths = generate_training_sets(PATH)
dst = '/Users/noa/Workspace/data/supermatrices_edited'
os.mkdir(dst)

for path in training_set_paths:
        try:
            with open(path) as original:
                reduced_local_alignment_data = list(SeqIO.parse(original, 'fasta'))
                n_seq = len(reduced_local_alignment_data)
                n_loci = len(reduced_local_alignment_data[0].seq)
                new_file_name = path.replace('/Users/noa/Workspace/data', "").replace('/', "_").replace('aln', 'fasta')

        except:
            try:
                with open(path) as original:
                    reduced_local_alignment_data = list(SeqIO.parse(original, 'phylip-relaxed'))
                    n_seq = len(reduced_local_alignment_data)
                    n_loci = len(reduced_local_alignment_data[0].seq)
                    new_file_name = path.replace('/Users/noa/Workspace/data', "").replace('/', "_").replace('aln',
                                                                                                            'phy')
            except:
                    print("got here")

        print(n_loci, n_seq)
        if n_loci >=10000 and n_seq>=15:
            shutil.copy2(path, os.path.join(dst,new_file_name))








