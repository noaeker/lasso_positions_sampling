
import pandas as pd
from help_functions import *
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def two_lines(x, a, b, c, d):
    one = a*x + b
    two = c*x+d
    return np.maximum(one, two)



lasso_path = ("/Users/noa/Workspace/lasso_positions_sampling_results/brlen_lasso_curr.csv")
output_data= pd.DataFrame(
)
unified_df = pd.read_csv(lasso_path)
curr_run_directory = "/Users/noa/Workspace/lasso_positions_sampling_results/test_dir"
create_dir_if_not_exists(curr_run_directory)
derivs = []
forward_derivs = []

for brlen_generator in unified_df["brlen_generator"].unique():
    print("brlen generator = {}".format(brlen_generator))
    res = unified_df.loc[(unified_df["brlen_generator"] == brlen_generator)]
    plateau = -1
    groupby_training_size = res.groupby("actucal_training_size").agg(lasso_test_R2=('lasso_test_R^2','mean')).reset_index()
    x = groupby_training_size["actucal_training_size"]
    y = groupby_training_size["lasso_test_R2"]
    x.diff = np.diff(x)
    y.diff = np.diff(y)
    deriv = y.diff / x.diff
    forward_deriv = [np.mean(deriv[i + 1:]) for i in range(len(deriv) - 1)]

    for i in range(len(deriv) - 1):
        if deriv[i] < forward_deriv[i] * 2 and deriv[i] > 0:
            plateau = x.iloc[i]
            break

    print(plateau)
    plt.scatter(x, y)
    plt.title(brlen_generator)
    plt.show()
    if output_data.empty:
        output_data = res
    else:
        output_data = pd.concat([output_data, res])








for dataset_id in unified_df["dataset_id"].unique():
    print("dataset id = {}".format(dataset_id))
    for brlen_generator in unified_df["brlen_generator"].unique():
        print("brlen generator = {}".format(brlen_generator))
        res = unified_df.loc[(unified_df["dataset_id"] == dataset_id) & (unified_df["brlen_generator"] == brlen_generator)]
        plateau=-1
        x = res["actucal_training_size"]
        y = res["lasso_test_R^2"]
        x.diff=np.diff(x)
        y.diff=np.diff(y)
        deriv = y.diff/x.diff
        forward_deriv = [np.mean(deriv[i+1:]) for i in range(len(deriv)-1)]

        for i in range(len(deriv)-1):
            if deriv[i]<forward_deriv[i]*2 and deriv[i]>0:
                plateau=x.iloc[i]
                break

        print(plateau)
        plt.scatter(x, y)
        plt.title(brlen_generator)
        plt.show()
        if output_data.empty:
            output_data= res
        else:
            output_data=pd.concat([output_data, res])
output_data.to_csv("test.csv")




output_data.to_csv("test.csv")

