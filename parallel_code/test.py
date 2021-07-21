import glmnet

import glmnet
import pandas as pd

X  = pd.read_csv("/Users/noa/Workspace/lasso_positions_sampling_results/training.csv")
y = X.sum(axis=1)
mfit = glmnet(x = X.copy(), y = y.copy(), family = 'mgaussian')
#fit = glmnet(x = X.copy(), y = y.copy(), family = 'gaussian', nlambda = 20)
