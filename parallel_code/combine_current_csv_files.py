
from help_functions import *

unified_csv_path = '/groups/pupko/noaeker/lasso_positions_sampling_results/current_spr.csv'
job_folder = '/groups/pupko/noaeker/lasso_positions_sampling_results/brlen_spr'


csv_paths= [os.path.join(job_folder,"job_{}".format(i),"{}.csv".format(i)) for i in range(1)]
all_MSA_results = add_csvs_content(csv_paths, unified_csv_path)
