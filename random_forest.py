from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
import pandas as pd
import numpy as np

X = pd.read_csv("samples.csv") # Input: parameter values
X = X.iloc[:, 1:] # Trim first column
y = pd.read_csv("all_tumour_counts.csv", header=None).to_numpy() # Output: tumour count
y = np.ravel(y) # Convert from column vector to 1D array

clf = RandomForestRegressor()
clf.fit(X,y)
importance = permutation_importance(clf, X, y)

sorted_indices = importance.importances_mean.argsort()

for i in reversed(sorted_indices):
    print(f"Parameter {i}: {importance.importances_mean[i]}")
