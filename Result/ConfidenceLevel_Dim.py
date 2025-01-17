import scipy.special as sp
from scipy.stats import norm

def sigma_to_p(n_sigma):
    """
    Converts the sigma level (nσ) to a probability (p).

    Parameters:
        n_sigma (float): The number of standard deviations (e.g., 1, 2, 3).

    Returns:
        float: The corresponding probability (p).
    """
    # Two-tailed probability for the given sigma
    p = 1 - 2 * norm.sf(n_sigma)
    return p

def find_delta(sigma=1, nu=1, tol=1e-10, max_iter=1e8):
    p = sigma_to_p(sigma)
    """
    Finds ∆χ² such that Q(ν/2, ∆/2) = 1 - p with high precision.

    Parameters:
        p (float): Desired probability.
        nu (float): Degrees of freedom.
        tol (float): Tolerance for convergence.
        max_iter (int): Maximum iterations for bisection.

    Returns:
        float: The value of ∆χ².
    """
    a = nu / 2.0  # Gamma function parameter
    target = 1.0 - p  # Target value for Q(a, ∆/2)
    
    # Initial bounds for bisection
    lower = 0.0
    upper = 10.0
    
    # Dynamically expand upper bound to ensure it brackets the root
    while sp.gammaincc(a, upper / 2.0) > target:
        upper *= 2.0  # Exponential growth for robust bracketing
    
    # Bisection method for root finding
    for _ in range(int(max_iter)):
        mid = (lower + upper) / 2.0
        value = sp.gammaincc(a, mid / 2.0)  # Regularized upper gamma Q(a, ∆/2)
        
        if abs(value - target) < tol:  # Convergence check
            return mid
        
        if value > target:  # Adjust bounds
            lower = mid
        else:
            upper = mid
    
    raise RuntimeError("Bisection method did not converge")

if __name__ == "__main__":
    # Define the desired sigma level and degrees of freedom
    sigma = 1
    nu = 4

    # Find the corresponding ∆χ² value
    delta = find_delta(sigma, nu)
    print(f"∆χ² for {sigma}σ and {nu} degrees of freedom: {delta}")