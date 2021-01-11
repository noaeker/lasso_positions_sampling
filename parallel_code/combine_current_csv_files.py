
from help_functions import *

unified_csv_path = '/groups/pupko/noaeker/lasso_positions_sampling_results/current_spr.csv'
job_folder = '/groups/pupko/noaeker/lasso_positions_sampling_results/brlen_spr'


csv_paths= [os.path.join(job_folder,"job_{}".format(i),"{}.csv") for i in range(50)]
all_MSA_results = add_csvs_content(csv_paths, unified_csv_path)
