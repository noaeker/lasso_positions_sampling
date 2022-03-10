
import re
import numpy as np
import os
import more_itertools as mit

def edit_num_locis_in_model_file_no_partition(model_path, n_loci):
    with open(model_path) as MODEL_FILE:
        text = MODEL_FILE.read()
    new_text = re.sub('noname = 1\-\d+',f'noname = 1-{n_loci}',text)
    with open(model_path,'w') as MODEL_FILE:
        MODEL_FILE.write(new_text)

def edit_frequency_synax_in_original_model_path(model_path):
    with open(model_path) as MODEL_FILE:
        text = MODEL_FILE.read()
    new_text = re.sub('[^+]F, ',f'+F, ',text)
    with open(model_path,'w') as MODEL_FILE:
        MODEL_FILE.write(new_text)


def generate_loci_corrected_partition_model_file(partition_results, partition_ind_to_name_dict, curr_run_directory, positions_subset = None):
    if positions_subset:
        partition_results = np.take(partition_results,positions_subset)
    model_file = os.path.join(curr_run_directory,'loci_adjusted_raxml_model')
    partition_indexes = np.unique(partition_results)
    with open(model_file,'w') as MODEL_FILE:
        for partition_ind in partition_indexes:
            partition_name = partition_ind_to_name_dict[partition_ind]
            MODEL_FILE.write(f"{partition_name} = ")
            partition_ind_indexes = list(np.where(partition_results == partition_ind)[0] + 1)
            ranges = [list(group) for group in mit.consecutive_groups(partition_ind_indexes)]
            ranges_str = ",".join([f"{range[0]}-{range[-1]}" for range in ranges])
            MODEL_FILE.write(f"{ranges_str}\n")
    return model_file



def parse_raxml_partition_file(model_file, orig_n_loci):
    per_loci_partition = np.zeros(orig_n_loci)
    partition_ind_to_name = {}
    with open(model_file) as MODEL_FILE:
        partitions = MODEL_FILE.readlines()
    for i,raw_partition in zip(range(1,len(partitions)+1),partitions):
        site_partitions  = re.findall('\s(\d+)\-(\d+)', raw_partition)
        partition_name = raw_partition.split("=")[0]
        partition_ind_to_name[i] = partition_name
        for site_group in site_partitions:
            start =int(site_group[0])-1
            end = int(site_group[1])
            np.put(per_loci_partition,ind =np.array(range(start,end)),v  = i)

    return per_loci_partition.astype(int), partition_ind_to_name
