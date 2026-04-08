import numpy as np

def Pmin(epsilon):
    return 1/(1+np.exp(-epsilon))
def Pmax(epsilon,Ko,Kc):
    return 1/(1+(Ko/Kc)**5*np.exp(-epsilon))
def PopenNorm(c,epsilon,Ko,Kc):
    return (Popen(c,epsilon,Ko,Kc) - Pmin(epsilon))/(Pmax(epsilon,Ko,Kc) - Pmin(epsilon))

def PmaxHetero(epsilon, Kos, Kcs):
    # Using list comprehension for element-wise division if inputs are lists/arrays
    ratio_prod = np.prod([Ko/Kc for Ko, Kc in zip(Kos, Kcs)])
    return 1/(1 + ratio_prod * np.exp(-epsilon))
def PoHetero(c, epsilon, Kos, Kcs):
    """
    Calculates absolute open probability for a heteropentamer.
    c: float or np.array of ligand concentrations
    epsilon: float, total gating energy
    Kos, Kcs: lists/arrays of length 5 containing dissociation constants
    """
    # Ensure c is handled correctly if it's a scalar or array
    term_o = np.prod([1 + c/Ko for Ko in Kos], axis=0)
    term_c = np.prod([1 + c/Kc for Kc in Kcs], axis=0)
    term_eps = np.sum(epsilon,axis=0)
    return term_o / (term_o + np.exp(-term_eps) * term_c)
def PoHeteroNorm_vectorized(c, epsilon, Kos, Kcs):
    """
    Vectorized calculation of the Normalized Open Probability for Hetero-pentamers.
    Uses the MWC state ratio form: P = 1 / (1 + L_iso * R)
    """
    # 1. Isomerization constant L (Closed/Open equilibrium at c=0)
    # epsilon shape: (L, Nr, 1, 1) -> L_iso: (L, Nr, 1)
    L_iso = np.exp(-epsilon).squeeze(-1)
    
    # 2. State Ratios R(c) = Prod_subunits [ (1+c/Kc) / (1+c/Ko) ]
    # Calculates the relative stability of Closed vs Open state for all c
    # Shape: (L, Nr, Nc)
    R_c = np.prod((1 + c / Kcs) / (1 + c / Kos), axis=2)
    
    # 3. Limiting Ratio R(inf) = Prod(Ko/Kc)
    # This defines the saturation limit as c -> infinity
    # Kos/Kcs shape: (L, Nr, 5, 1) -> prod axis 2 -> (L, Nr, 1)
    # We DO NOT squeeze the last dim so it broadcasts with L_iso (L, Nr, 1)
    R_inf = np.prod(Kos / Kcs, axis=2)
    
    # 4. Calculate Probabilities: P = 1 / (1 + L * R)
    P_c = 1.0 / (1.0 + L_iso * R_c)
    P_min = 1.0 / (1.0 + L_iso)       # At c=0, R=1, shape (L, Nr, 1)
    P_max = 1.0 / (1.0 + L_iso * R_inf) # Shape (L, Nr, 1)
    
    # 5. Normalize
    # Use np.divide with 'where' to handle cases where P_max == P_min (flat response)
    delta_P = P_max - P_min
    return np.divide(P_c - P_min, delta_P, where=(delta_P != 0), out=np.zeros_like(P_c))

def PoHeteroNorm(c, epsilon, Kos, Kcs):
    p_abs = PoHetero(c, epsilon, Kos, Kcs)
    p_min = Pmin(epsilon)
    p_max = PmaxHetero(epsilon, Kos, Kcs)
    # Normalize
    return (p_abs - p_min) / (p_max - p_min)

def generate_discrete_curve(epsilon, Kos, Kcs, Nc, c_min, c_max, Na):
    """
    Generates a discrete activation curve from the continuous Popen model.

    Args:
        epsilon, Ko, Kc: Parameters for the Popen function.
        Nc (int): Number of concentration bins.
        c_min (float): Minimum concentration for binning.
        c_max (float): Maximum concentration for binning.
        Na (int): Number of activity bins.

    Returns:
        np.array: A 1D array of length Nc containing the discrete
                  activity curve (e.g., [0, 0, 1] for Na=2, Nc=3).
    """
    
    # --- Concentration Bins (Logarithmic) ---
    # Create Nc+1 bin edges from c_min to c_max
    c_edges = np.logspace(np.log10(c_min), np.log10(c_max), Nc + 1)
    
    # Get Nc bin representative centers (geometric mean)
    c_centers = np.sqrt(c_edges[:-1] * c_edges[1:])
    
    # --- Activity Bins (Linear) ---
    # Create Na+1 bin edges from 0 to 1
    a_edges = np.linspace(0, 1, Na + 1)
    a_centers = 0.5*(a_edges[1:] + a_edges[:-1])
    a_centers[0] = 0.
    a_centers[-1] = 1.
    
    # --- Discretization ---
    # 1. Get the continuous Popen value at each concentration center
    #continuous_values = PopenNorm(c_centers, epsilon, Ko, Kc)
    continuous_values = PoHeteroNorm(c_centers, epsilon, Kos, Kcs)
    
    # 2. Digitize these values into bins
    discrete_curve = np.digitize(continuous_values, a_edges)
    return c_edges,a_centers[discrete_curve-1]

def generate_continuous_A(epsilon,Kos,Kcs,c_centers,a_edges,indices_to_avg=None):
    # --- 1. Prepare Input Shapes for Broadcasting ---
    
    # c: Add axes for L, Nr, and Subunit -> (1, 1, 1, Nc)
    c_grid = c_centers[np.newaxis, np.newaxis, np.newaxis, :]
    
    # Kos, Kcs: Add axis for Concentration -> (L, Nr, 5, 1)
    Kos_grid = Kos[:, :, :, np.newaxis]
    Kcs_grid = Kcs[:, :, :, np.newaxis]
    
    # epsilon: All 5 subunits share the same epsilon for the complex global transition.
    # We take the first value and add axes -> (L, Nr, 1, 1)
    eps_grid = epsilon[:, :, 0][:, :, np.newaxis, np.newaxis]
    
    # --- 2. Calculate Continuous Response ---

    continuous_A = PoHeteroNorm_vectorized(c_grid, eps_grid, Kos_grid, Kcs_grid)

    if indices_to_avg is not None:
        all_indices = np.arange(continuous_A.shape[1])
        remaining_indices = np.delete(all_indices, indices_to_avg) # [0, 1]
        # 2. Extract the "untouched" rows
        untouched_part = continuous_A[:, remaining_indices, :]
        # 3. Calculate the mean of the "target" rows
        # We use axis=1 to average across the row dimension
        avg_part = np.mean(continuous_A[:, indices_to_avg, :], axis=1, keepdims=True)
        # 4. Combine them
        continuous_A = np.concatenate([untouched_part, avg_part], axis=1)

    return continuous_A
    

def generate_activity_matrix(epsilon, Kos, Kcs, c_centers, a_edges,indices_to_avg=None):
    """
    Generates the discrete activity matrix for a set of ligands.
    
    Args:
        epsilon: (L, Nr, 5)
        Kos:     (L, Nr, 5)
        Kcs:     (L, Nr, 5)
        c_centers: (Nc,)
        a_edges:   (Na+1,)
        
    Returns:
        discrete_A: Integer array of shape (L, Nr, Nc)
    """
    L, Nr, penta = epsilon.shape
    Nc = len(c_centers)
    
    # --- 2. Calculate Continuous Response ---
    continuous_A = generate_continuous_A(epsilon,Kos,Kcs,c_centers,a_edges,indices_to_avg)
    
    # --- 3. Discretize ---
    # Clip to ensure numerical stability doesn't push values < 0 or > 1
    continuous_A = np.clip(continuous_A, 0, 1)

    # Digitize using inner edges as thresholds
    discrete_A = np.digitize(continuous_A, a_edges[1:-1])
    
    return discrete_A

import itertools
import random

def generate_Ks(L=100, n_types=2, n_subunits=5, n_hetero_sample='all', Kc_Amp=10., eps_Amp=-10.,K_generation = 'uniform'):
    """
    Generates Ks, Kos, and Epsilon for N types with memory-efficient heteromer selection.
    
    Args:
        L (int): Number of simulations.
        n_types (int): Number of base subunit types.
        n_subunits (int): Number of subunits per receptor.
        n_hetero_sample (int, str, list, or dict): 
            - 'all': Keep all heteromers.
            - int (e.g., 100): Randomly select this many heteromers (any size).
            - dict (e.g., {2: 100, 3: 50}): Select 100 combos of size 2, 50 of size 3, etc.
            - list (e.g., [0, 0, 100, 50]): Index i corresponds to set size i.
        K_generation (str) : how K values are generated, notice that for a good ion channel
                            we need to keep forall ligand and homomer : Ko < Kc : dissociation
                            smaller when open that when closed            
            - 'uniform' each Kc and Ko are drawn from a uniform distrib
            - 'hierarchical' the first Kc, Ko are drawn from uniform distrib, next one is smaller/bigger etc...
    """
    
    # 1. Generate Homomers directly (Always keep all homomers)
    # This is fast and tiny in memory
    homomers = [(i,) * n_subunits for i in range(n_types)]
    
    # 2. Setup Reservoir Sampling for Heteromers
    # This determines how we store data while iterating
    mode = 'all'
    reservoirs = {} # Stores the selected combos
    counts = {}     # Tracks how many we've seen (for probability calculation)
    
    if n_hetero_sample == 'all':
        mode = 'all'
        reservoirs['all'] = []
    elif isinstance(n_hetero_sample, int):
        mode = 'total'
        reservoirs['total'] = []
        counts['total'] = 0
        limit_total = n_hetero_sample
    elif isinstance(n_hetero_sample, (list, dict)):
        mode = 'stratified'
        # Normalize input to a dict: {size: limit}
        if isinstance(n_hetero_sample, list):
            # mapping index -> limit
            limits = {i: val for i, val in enumerate(n_hetero_sample) if val > 0}
        else:
            limits = n_hetero_sample
            
        # Initialize reservoirs for each requested size
        for size in limits:
            reservoirs[size] = []
            counts[size] = 0
            
    # 3. Iterate via Generator (Memory Efficient)
    # We DO NOT convert this to a list. We process one by one.
    all_combos_iter = itertools.combinations_with_replacement(range(n_types), n_subunits)
    
    for c in all_combos_iter:
        unique_count = len(set(c))
        
        # Skip homomers (handled separately)
        if unique_count == 1:
            continue
            
        # --- Mode: Keep All ---
        if mode == 'all':
            reservoirs['all'].append(c)
            
        # --- Mode: Total Random Limit ---
        elif mode == 'total':
            # Standard Reservoir Sampling
            counts['total'] += 1
            current_count = counts['total']
            
            if len(reservoirs['total']) < limit_total:
                reservoirs['total'].append(c)
            else:
                # Randomly replace an existing item with diminishing probability
                r = random.randint(0, current_count - 1)
                if r < limit_total:
                    reservoirs['total'][r] = c
                    
        # --- Mode: Stratified (Specific counts for specific sizes) ---
        elif mode == 'stratified':
            if unique_count in limits:
                counts[unique_count] += 1
                current_count = counts[unique_count]
                limit = limits[unique_count]
                
                if len(reservoirs[unique_count]) < limit:
                    reservoirs[unique_count].append(c)
                else:
                    r = random.randint(0, current_count - 1)
                    if r < limit:
                        reservoirs[unique_count][r] = c

# 4. Flatten and Sort Results by Complexity (Size of Set)
    selected_hetero = []

    if mode == 'all' or mode == 'total':
        # Combine lists if needed (usually just one exists)
        raw_heteros = reservoirs.get('all', []) + reservoirs.get('total', [])
        
        # Sort first by "number of unique types" (complexity), then numerically
        selected_hetero = sorted(raw_heteros, key=lambda x: (len(set(x)), x))

    elif mode == 'stratified':
        # Loop through keys in order (2, 3, 4...) to maintain complexity order
        for size in sorted(reservoirs.keys()):
            # Sort the combos within that specific size bucket numerically
            bucket_combos = sorted(reservoirs[size])
            selected_hetero.extend(bucket_combos)

    # Homomers are already sorted (0,0..), (1,1..)
    homomers.sort()
    
    # Final Concatenation: Homomers -> Hetero(2) -> Hetero(3) -> ...
    final_combos = homomers + selected_hetero
    
    Nr = len(final_combos)
    
    # 5. Pre-allocate arrays
    Kcs = np.zeros((L, Nr, n_subunits), dtype=float)
    Kos = np.zeros((L, Nr, n_subunits), dtype=float)
    epsilon = np.zeros((L, Nr, n_subunits), dtype=float)
    
    # 6. Generate Parameters
    for l in range(L):
        # --- A. Generate Base Parameters for N types ---
        
        if K_generation=='hierarchical':
            base_Kc,base_Ko,base_eps = [],[],[]
            prev_Kc = 0
            prev_Ko = None        
            for i in range(n_types):
                # Kc increases
                current_Kc = (prev_Kc if i > 0 else 0) + np.random.random() * Kc_Amp
                base_Kc.append(current_Kc)
                
                # Ko decreases (and Ko <= Kc)
                if i == 0:
                    current_Ko = np.random.random() * current_Kc
                else:
                    current_Ko = np.random.random() * prev_Ko
                base_Ko.append(current_Ko)
                
                base_eps.append(np.random.random() * eps_Amp)
                
                prev_Kc = current_Kc
                prev_Ko = current_Ko
        elif K_generation == 'uniform':
            base_Kc = [np.random.random() * Kc_Amp for _ in range(n_types)]
            base_Ko = [np.random.random() * Kc for Kc in base_Kc]
            base_eps = [np.random.random() * eps_Amp for _ in range(n_types)]

        # --- B. Map to Receptor Structure ---
        for r_idx, combo_indices in enumerate(final_combos):
            Kcs[l, r_idx, :] = [base_Kc[i] for i in combo_indices]
            Kos[l, r_idx, :] = [base_Ko[i] for i in combo_indices]
            
            combo_eps_values = [base_eps[i] for i in combo_indices]
            epsilon[l, r_idx, :] = [np.mean(combo_eps_values)] * n_subunits
            
    return Kcs, Kos, epsilon
