import numpy as np

def ec50_hetero_from_homo(c12_homo1,c12_homo2,k):
    return c12_homo1**k * c12_homo2**(5-k)


def ec50_approx(epsilon, Kos):
    # Geometric mean of Ko * exp(-epsilon/5)
    Ko_geom = np.prod(Kos)**(1/5)
    epsilon_tot = np.sum(epsilon)
    return Ko_geom * np.exp(-epsilon_tot/5)


def get_exact_ec50(epsilon, Kos, Kcs):
    """
    Finds the concentration c where PoHeteroNorm(c) = 0.5
    Uses root finding for precision.
    """
    # Target absolute probability corresponding to Normalized=0.5
    p_min = Pmin(epsilon)
    p_max = PmaxHetero(epsilon, Kos, Kcs)
    target_p = 0.5 * (p_max - p_min) + p_min
    
    # Define function to minimize: P_calc - P_target
    def func(log_c):
        c = 10**log_c
        return PoHetero(c, epsilon, Kos, Kcs) - target_p
    
    # Search in a wide log range (e.g., 1e-12 to 1e0)
    # We create a dynamic bracket based on the approximation to ensure convergence
    guess = EC50_approx(epsilon, Kos)
    log_guess = np.log10(guess)
    
    try:
        root_log = brentq(func, log_guess - 5, log_guess + 5)
        return 10**root_log
    except ValueError:
        return np.nan # Solver failed (likely flat curve or extreme params)