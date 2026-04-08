import numpy as np
def response_to_index(response_vector, Na):
    """
    Maps a response vector [a0, a1, ..., aNr-1] to a unique integer index.
    This treats the vector as a number in base Na.
    """
    index = 0
    power = 1
    for activity_level in response_vector:
        index += activity_level * power
        power *= Na
    return int(index)
def response_to_index_vectorized(response_matrix, Na):
    """
    Maps a matrix of response vectors to unique integer indices.
    
    Args:
        response_matrix: (N_events, Nr) array. Each row is one event [r0, r1, r2].
        Na: Number of activity bins (base).
    Returns:
        indices: (N_events,) array of unique integers.
    """
    Nr = np.asarray(response_matrix).shape[1]
    # Powers array: [1, Na, Na^2, ...]
    powers = Na ** np.arange(Nr)
    # Dot product broadcasting to get unique index for each row
    return np.sum(response_matrix * powers, axis=1)
def index_to_response(index, Nr, Na):
    """Decodes unique index back to vector [r0, r1, r2...]"""
    response = np.zeros(Nr, dtype=int)
    current = index
    for i in range(Nr):
        response[i] = current % Na
        current //= Na
    return response

# ----------------------------------------------------------------
# -----------Uncomment to perform tests on the indexing-----------
# ----------------------------------------------------------------
# Write the response of three receptors. They must be inside activity bins
#y_index = [[0,Na//2, Na-1]]
#y = [a_centers[i] for i in y_index]
#print([np.digitize(y_val,a_centers)-1 for y_val in y])
#print(y)
#print(response_to_index(y_index[0],Na))
#print(response_to_index_vectorized(y_index,Na))
#print(index_to_response(response_to_index(y_index[0],Na),Nr=3,Na=Na))

from scipy.special import entr

def compute_mi(discrete_A):

    L, Nr, Nc = discrete_A.shape
    
    N_total = L * Nc
    
    events_matrix = discrete_A.transpose(0, 2, 1).reshape(-1, Nr)
    
    _, counts = np.unique(events_matrix, axis=0, return_counts=True)
    
    # M_A contains the counts of only the states that actually occurred.
    M_A = counts
    
    # 4. Compute Probabilities P(A)
    P_A = M_A / N_total
    
    MI = np.sum(entr(P_A))
    
    return MI
def compute_P(discrete_A,Na):
    L, Nr, Nc = discrete_A.shape
    
    N_total = L * Nc
    
    events_matrix = discrete_A.transpose(0, 2, 1).reshape(-1, Nr)
    
    values, counts = np.unique(events_matrix, axis=0, return_counts=True)
    
    Pa = np.zeros(([Na for _ in range(Nr)]),dtype=float)

    for value,count in zip(values,counts):
        Pa[tuple(value)] = count
    Pa /= N_total

    return Pa

def compute_P_sparse(discrete_A):
    L, Nr, Nc = discrete_A.shape
    N_total = L * Nc
    
    events_matrix = discrete_A.transpose(0, 2, 1).reshape(-1, Nr)
    
    # Get unique rows and their counts
    values, counts = np.unique(events_matrix, axis=0, return_counts=True)
    
    # Calculate probabilities for observed states only
    probs = counts / N_total
    
    # Return a dictionary mapping {State_Tuple: Probability}
    # distinct_states is a list of tuples, e.g., [(0,1,0...), (1,1,0...)]
    distinct_states = [tuple(v) for v in values]
    
    return dict(zip(distinct_states, probs))


# --- Helper to Compute Marginal P(a) ---
def compute_marginal_P(discrete_A, Na):
    """
    Computes the marginal probability of each activity bin (0 to Na-1)
    observed across all ligands and receptors in the matrix.
    
    Returns:
        prob_vector: Array of shape (Na,) summing to 1.
    """
    # Flatten to treat all observations equally
    flat_A = discrete_A.ravel()
    # Count occurrences of each bin
    counts = np.bincount(flat_A, minlength=Na)
    # Normalize
    return counts / flat_A.size
