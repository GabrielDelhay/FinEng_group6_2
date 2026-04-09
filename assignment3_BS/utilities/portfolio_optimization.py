"""
Portfolio optimization utilities.

Implements:
- Minimum variance portfolio (closed-form solution)
- Mean-variance portfolio (closed-form solution)

References:
- Markowitz, H. (1952). "Portfolio Selection." The Journal of Finance.
"""

import numpy as np
from utilities.covariance_utilities import (
    _validate_covariance_matrix,
)


def minimum_variance_portfolio(cov_matrix: np.ndarray) -> np.ndarray:
    """
    Calculate the minimum variance portfolio weights given a covariance matrix.
    In particular the weights are given by: (pay attention: the following formula has to be written with the invese ov matrix)
    w = (Σ * 1) / (1^T * Σ * 1), i.e. the solution of the optimization problem: 
    min_w w^T * Σ * w, subject to 1^T * w = 1.

    Parameters:
        cov_matrix (np.ndarray): Covariance matrix of asset returns.

    Returns:
        np.ndarray: Weights of the minimum variance portfolio.
    """
    cov_matrix = _validate_covariance_matrix(
        cov_matrix,
        name="cov_matrix",
        require_positive_definite=True,
        positive_definite_message=(
            "cov_matrix must be positive definite (symmetric with positive eigenvalues)"
        ),
    )

    n = cov_matrix.shape[0]
    ones_vec = np.ones((n, 1))

    min_var_ptf_numerator = np.linalg.inv(cov_matrix) @ ones_vec
    min_var_ptf_weights = min_var_ptf_numerator / (ones_vec.T @ np.linalg.inv(cov_matrix) @ ones_vec)

    return min_var_ptf_weights.flatten()


def mean_variance_portfolio(
    expected_returns: np.ndarray,
    cov_matrix: np.ndarray,
    risk_aversion: float = 1.0,
) -> np.ndarray:
    """
    Calculate the classic mean-variance portfolio weights given expected returns and a
    covariance matrix.

    In particular the weights solve:
    max_w mu^T * w - (gamma / 2) * w^T * Sigma * w, subject to 1^T * w = 1,
    where mu are the expected returns and gamma is the risk-aversion parameter.

    Parameters:
        expected_returns (np.ndarray): Expected returns vector.
        cov_matrix (np.ndarray): Covariance matrix of asset returns.
        risk_aversion (float): Risk-aversion parameter gamma. Must be strictly positive.

    Returns:
        np.ndarray: Weights of the mean-variance portfolio.
    """
    cov_matrix = _validate_covariance_matrix(
        cov_matrix,
        name="cov_matrix",
        require_positive_definite=True,
        positive_definite_message=(
            "cov_matrix must be positive definite (symmetric with positive eigenvalues)"
        ),
    )

    expected_returns = np.asarray(expected_returns, dtype=float)
    if expected_returns.ndim == 2 and 1 in expected_returns.shape:
        expected_returns = expected_returns.reshape(-1)
    elif expected_returns.ndim != 1:
        raise ValueError(
            "expected_returns must be one-dimensional or a single-column vector"
        )

    if expected_returns.shape[0] != cov_matrix.shape[0]:
        raise ValueError(
            "expected_returns and cov_matrix must refer to the same number of assets, "
            f"got {expected_returns.shape[0]} and {cov_matrix.shape[0]}"
        )

    if not np.isfinite(expected_returns).all():
        raise ValueError("expected_returns contains NaN or Inf values")

    if not np.isfinite(risk_aversion):
        raise ValueError("risk_aversion must be finite")

    if risk_aversion <= 0:
        raise ValueError(
            f"risk_aversion must be strictly positive, got {risk_aversion}"
        )
    n = cov_matrix.shape[0]
    ones_vec = np.ones(n)

    # Σ⁻¹μ e Σ⁻¹1
    inv_cov_mu = np.linalg.solve(cov_matrix, expected_returns)
    inv_cov_ones = np.linalg.solve(cov_matrix, ones_vec)

    A = ones_vec @ inv_cov_mu      # 1ᵀ Σ⁻¹ μ
    C = ones_vec @ inv_cov_ones    # 1ᵀ Σ⁻¹ 1

    # a* = (A/γ) · (Σ⁻¹μ / A) + (1 - A/γ) · (Σ⁻¹1 / C)
    mean_var_ptf_weights = (A / risk_aversion) * (inv_cov_mu / A) + (1 - A / risk_aversion) * (inv_cov_ones / C)
    return mean_var_ptf_weights.flatten()

# add a bisection funciotn in order to get to the right gamma
def bisection(f, a, b, tol=1e-6, max_iter=100):
    """
    Finds a root of the function f in the interval [a, b] using the bisection method.
    
    Parameters:
    f (callable): The function for which to find the root.
    a (float): Left endpoint of the interval.
    b (float): Right endpoint of the interval.
    tol (float): Tolerance for the stopping criterion (interval width).
    max_iter (int): Maximum number of allowed iterations.
    
    Returns:
    float: An approximation of the root.
    """
    # Initial check: f(a) and f(b) must have opposite signs
    if f(a) * f(b) >= 0:
        raise ValueError("The function must have opposite signs at the endpoints of the interval [a, b].")
    
    for i in range(max_iter):
        # Calculate the midpoint
        c = (a + b) / 2.0
        
        # Stopping criterion: 
        # Exact root found or the interval has become smaller than the tolerance
        if f(c) == 0 or (b - a) / 2.0 < tol:
            return c
        
        # Narrow down the interval for the next iteration
        if f(a) * f(c) < 0:
            b = c  # The root is between 'a' and 'c'
        else:
            a = c  # The root is between 'c' and 'b'
            
    print("Warning: maximum number of iterations reached without satisfying the tolerance.")
    return (a + b) / 2.0
