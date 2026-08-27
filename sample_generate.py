import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Settings
parameters_file = "parameters_sensitivity_analysis.csv"
size = 200 # Total number of simulation runs - change to 5000
num_par = 25 # Number of parameters
samples_df = pd.DataFrame()

# Uniform variation of parameters
'''
# Read parameters file
with open(parameters_file, newline='') as f:
    reader = csv.reader(f)

    next(reader) # Skip the header

    # Generate samples according to a normal distribution
    for i, line in enumerate(reader):
        
        mean = float(line[1])
        samples = np.random.default_rng().uniform(low=0.25*mean, high=1.75*mean, size=size)

        samples_column = np.array(samples)
        samples_df[f"parameter_{i}"] = samples_column # Column of df is titled as "parameter #"
'''

# Lognormal variation of parameters

# Convert from normal to lognormal means and sds
def norm_to_lognorm(mean, sd):
    var = sd**2

    mu = np.log(mean / (np.sqrt(1+(var/(mean**2)))))
    sigma = np.sqrt(np.log(1+var/(mean**2)))

    return mu, sigma

# Read parameters file
with open(parameters_file, newline='') as f:
    reader = csv.reader(f)

    next(reader) # Skip the header

    # Generate samples according to a normal distribution
    for i, line in enumerate(reader):

        mean = float(line[1]) # Mean of the underlying normal distribution
        sd = float(line[2]) # SD of the underlying normal distribution
        low = float(line[3]) # Lower bound
        high = float(line[4]) # Upper bound

        print(mean)

        mu, sigma = norm_to_lognorm(mean, sd)

        samples = np.random.default_rng().lognormal(mean=mu, sigma=sigma, size=size)
        
        samples_column = np.array(samples)
        samples_df[f"parameter_{i}"] = samples_column # Column of df is titled as "parameter #"

        
samples_df.to_csv("samples.csv")