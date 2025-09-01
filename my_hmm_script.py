import numpy as np
from hmmlearn import hmm

# Define the HMM states
states = ["X", "Y", "Z"]
n_states = len(states)

# Transition probability matrix
t = 0.01
w = 0.0001

transition_matrix = np.array([
    [1 - t, t, 0],   # X -> X, X -> Y, X -> Z
    [w, 1 - 2*w, w], # Y -> X, Y -> Y, Y -> Z
    [0, t, 1 - t]    # Z -> X, Z -> Y, Z -> Z
])

# Mean and variance for emission distributions (Gaussian)
means = np.array([[-2.0], [0.0], [2.0]])  # Means for X, Y, Z
covars = np.array([[1.0], [1.0], [1.0]])  # Variance (std dev squared)

# Initialize the HMM model
model = hmm.GaussianHMM(n_components=n_states, covariance_type="diag", n_iter=1000)
model.startprob_ = np.array([1/3, 1/3, 1/3])  # Equal probability to start in any state
model.transmat_ = transition_matrix
model.means_ = means
model.covars_ = covars

# Observation sequence
observations = np.array([[-3], [1], [-3], [-2], [5], [1], [2], [1], [6], [1], [5], [4]])

# Run Viterbi Algorithm to find the most likely state sequence
logprob, hidden_states = model.decode(observations, algorithm="viterbi")

# Convert state indices to state names
decoded_path = [states[state] for state in hidden_states]

# Display the results
import pandas as pd

df = pd.DataFrame({
    "Observation": observations.flatten(),
    "Hidden State": decoded_path
})

print(df)
