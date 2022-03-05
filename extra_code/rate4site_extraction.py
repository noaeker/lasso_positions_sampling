from rate4site import *


def generate_commands(name, job_id, suffix):
    full_r4s_path = (
        f"/groups/pupko/noaeker/lasso_positions_sampling_results/f_lasso30/job_{job_id}/_groups_pupko_noaeker_data_supermatrices_edited__supermatrices_{name}_{name}.{suffix}/n_seq_30/n_loci_80000/r4s.res")
    local_full_r4s = f"/Users/noa/Workspace/lasso_positions_sampling_results/{name}_full_rate4site"
    command_full = (f"scp noaeker@power9login.tau.ac.il:{full_r4s_path} {local_full_r4s} ")
    #os.rename(full_r4s_path, f"{full_r4s_path}_{name}")
    lasso_r4s_path = f"/groups/pupko/noaeker/lasso_positions_sampling_results/f2_mg30/job_{job_id}/_groups_pupko_noaeker_data_supermatrices_edited__supermatrices_{name}_{name}.{suffix}/n_seq_30/n_loci_80000/Lasso_folder/exponential/training_4000_random_tree_eval/trimmed_4000/lasso_rate_4_site "
    #os.rename(lasso_r4s_path, f"{lasso_r4s_path}_{name}")
    local_lasso_r4s = f"/Users/noa/Workspace/lasso_positions_sampling_results/{name}_lasso_rate4site"

    command_lasso =f"scp noaeker@power9login.tau.ac.il:{lasso_r4s_path} {local_lasso_r4s}  "
    print(command_full)
    print("\n")
    print(command_lasso)
    print("\n")
    return local_full_r4s,local_lasso_r4s

def analayze_folder_rate4site(full_r4s_local_path,lasso_r4s_local_path):

    full_rate4site_numbers = parse_rate4site(full_r4s_local_path)
    with open(full_r4s_local_path,'w') as RATE4SITE:
        RATE4SITE.write("Numbers\n")
        RATE4SITE.writelines([str(n)+"\n" for n in full_rate4site_numbers])

    with open(lasso_r4s_local_path) as LASSO_RATE4SITE:
        lasso_rate4site_numbers = LASSO_RATE4SITE.read().split(' ')
    with open(lasso_r4s_local_path,'w') as LASSO_RATE4SITE:
         LASSO_RATE4SITE.write("Numbers\n")
         LASSO_RATE4SITE.writelines([str(n)+"\n" for n in lasso_rate4site_numbers])

def main():
        #for name, job_ind,suffix in zip(("NagyA1","ShenA9","YangA8","MisoA2","StruA5","WickA3"),(0,0,1,1,2,2), ('fasta','fasta','fasta','fasta','fasta','phy')):
        for name, job_ind, suffix in zip(("NagyA1",),
                                         (0,), ('fasta',)):
            full_r4s_local_path, lasso_r4s_local_path = generate_commands(name, job_ind, suffix)
            print(full_r4s_local_path, lasso_r4s_local_path)
            analayze_folder_rate4site(full_r4s_local_path, lasso_r4s_local_path)


if __name__ == "__main__":
    main()
