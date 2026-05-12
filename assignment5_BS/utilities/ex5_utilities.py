
import numpy as np
import pandas as pd
import datetime as dt
from typing import Tuple, List
from scipy.integrate import quad    # compute numerical integration

from utilities.ex1_utilities import (
    get_discount_factor_by_zero_rates_linear_interp,
    year_frac_act_x,
)

def simulate_ou_process(
    simulation_grid: pd.DatetimeIndex,
    today: pd.Timestamp,
    mean_reversion: float,
    sigma: float,
    n_sim: int,
    x0=0.0
) -> np.ndarray :
    """
    Simulate the Ornstein–Uhlenbeck (Hull–White) state process x_t.

    Parameters
    ----------
    simulation_grid : pd.DatetimeIndex
        Simulation dates
    today : pd.Timestamp
        Initial date
    mean_reversion : float
        Mean reversion parameter a
    sigma : float
        Volatility parameter
    n_sim : int
        Number of Monte Carlo scenarios
    x0 : float, default 0.0
        Initial value of the process

    Returns
    -------
    np.ndarray
        Simulated paths of x_t, shape (n_sim, n_dates)
    """

    # !!! COMPLETE AS APPROPRIATE !!! 
    
    n_steps = None

    # Generate standard normal variables
    Z = None

    # Allocate paths
    x_paths = np.zeros((n_sim, n_steps))
    x_prev = np.full(n_sim, x0)

    # Time stepping
    for i, sim_date in enumerate(simulation_grid):

        prev_date = None
        dt = year_frac_act_x(prev_date, sim_date, 360)

        x_current = None

        x_paths[:, i] = x_current
        x_prev = x_current

    return x_paths


def affine_trick(
    valuation_date: dt.datetime,
    pricing_grid: List[dt.datetime],
    mean_reversion: float,
    sigma: float,
    discount_factors: pd.Series,
) -> Tuple[pd.Series, pd.Series]:
    """
    Affine trick: Exploits the affine structure of the Hull-White model. Output the functions
    A(t, t_j) and C(t, t_j) pre-computed on the pricing grid s.t.
    B(t, t_j) = A(t, t_j) * exp(-C(t, t_j) * x(t)).

    Parameters:
        t (dt.datetime): Date w.r.t. which the computations are made.
        pricing_grid (List[dt.datetime]): Pricing grid.
        mean_reversion (float): Hull-white mean reversion speed.
        sigma (float): Hull-white interest rate volatility.
        discount_factors (pd.Series): Discount factors curve.

    Returns:
        Tuple[pd.Series, pd.Series]: Tuple with the precomputed functions A(t, t_j) and C(t, t_j).
    """
    
    # Reference date and containers
    reference_date = discount_factors.index[0].to_pydatetime()

    A = pd.Series(index=pricing_grid, dtype=float)
    C = pd.Series(index=pricing_grid, dtype=float)

    # Ensure discount factors exist on the whole pricing grid
    missing_dates = [d for d in pricing_grid if d not in discount_factors.index]
  
    # !!! COMPLETE AS APPROPRIATE !!! 
    if missing_dates:
        interpolated_dfs = None
        
        discount_factors = (
            pd.concat([discount_factors, interpolated_dfs])
            .sort_index()
        )

    # Volatility integrand
    sigma_func = lambda u: None

    # Compute A(t,T) and C(t,T)
    for T in pricing_grid:

        if T <= valuation_date:
            A[T] = 1
            C[T] = 0
            continue

        tau_t_T = None
        tau_0_T = None
        tau_0_t = None

        A[T] = None
        C[T] = None

    return A, C


    
def update_collateral(
    mtm: np.ndarray,
    VM: np.ndarray,
    THR_B: float = 0.0,
    THR_C: float = 0.0,
    MTA_B: float = 0.0,
    MTA_C: float = 0.0,
    Cap_B: float = np.inf,
    Cap_C: float = np.inf,
) -> np.ndarray:
    """
    Update the Variation Margin (VM) at a single time step under a bilateral CSA,
    from the perspective of party B, for multiple Monte Carlo scenarios.

    The function applies:
    - computation of target variation margin VM*,
    - margin calls subject to minimum transfer amounts (MTA),
    - application of collateral caps.

    Independent amount (IA), haircuts, rounding rules, settlement lags,
    collateral revaluation and remuneration effects are ignored here.

    Parameters
    ----------
    mtm : np.ndarray
        Array of mark-to-market values for each Monte Carlo scenario,
        measured from the perspective of party B.
        Shape: (n_scenarios,)

    VM : np.ndarray
        Array of variation margin currently posted for each scenario.
        Shape: (n_scenarios,)

    THR_B : float, default 0.0
        Threshold granted to party B (applies when MtM < 0).

    THR_C : float, default 0.0
        Threshold granted to counterparty C (applies when MtM > 0).

    MTA_B : float, default 0.0
        Minimum Transfer Amount when B is required to post collateral.

    MTA_C : float, default 0.0
        Minimum Transfer Amount when C is required to post collateral.

    Cap_B : float, default +inf
        Maximum collateral that B can be required to post.

    Cap_C : float, default +inf
        Maximum collateral that C can be required to post.

    Returns
    -------
    np.ndarray
        Updated variation margin for each scenario.
        Shape: (n_scenarios,)

    Notes
    -----
    Target variation margin is defined scenario-wise as:

        VM* =
            max(MtM - THR_C, 0)    if MtM >= 0
            min(MtM + THR_B, 0)    if MtM < 0

    Margin calls are applied only when the required change exceeds
    the relevant MTA.
    """

    # Compute target VM* (before any constrain)

    VM_star = None
    
    # Theoretical margin call

    delta_VM = None

    # Effective margin call (compare with MTA)
    
    margin_call = None
    
    # Update VM before caps

    VM_no_cap = None
    
    # Apply collateral caps

    VM_updated = None

    return VM_updated