
import pandas as pd
import logging
from sklearn import linear_model

training_set_path="/Users/noa/Workspace/lasso_positions_sampling_results/training_3200_exponential.csv"
test_set_path = "/Users/noa/Workspace/lasso_positions_sampling_results/test_100_exponential.csv"
training_set = pd.read_csv(training_set_path)
test_set= pd.read_csv(test_set_path)
test_r_squared= []
for n in [100,200,400,800,1600,3200]:
    for i in range(3):
        samp = training_set.sample(n=n, replace=False, random_state=1)
        y_training = samp.sum(axis=1)
        print("Sitelh df dimensions, for lasso computations, are: " + str(samp.shape))
        lasso_model = linear_model.LassoCV(cv=5, normalize=True, max_iter=100000, positive=True).fit(samp,
                                                                                                     y_training)  # add positive=True if using RaxML
        y_test = test_set.sum(axis=1)
        test_predicted_values = lasso_model.predict(test_set)
        test_r_squared = lasso_model.score(test_set, y_test)
        print("test_r_squared={}".format(test_r_squared))