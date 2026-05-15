
import numpy as np
import pandas as pd
import datetime as dt
from typing import Tuple, List
from scipy.integrate import quad    # compute numerical integration

from utilities.ex0_utilities import (                         # MISTAKE: it was imported from the wrong utilities 
    get_discount_factor_by_zero_rates_linear_interp,
    year_frac_act_x,
    year_frac_30e_360,
)

from utilities.ex1_utilities import (
    swap_par_rate,
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

    n_steps = len(simulation_grid)
    
    # Generate standard normal variables
    Z = np.random.standard_normal((n_sim, n_steps))

    # Allocate paths
    x_paths = np.zeros((n_sim, n_steps))
    x_prev = np.full(n_sim, x0)

    # Time stepping
    for i, sim_date in enumerate(simulation_grid):

        prev_date = today if i == 0 else simulation_grid[i - 1]
        dt = year_frac_act_x(prev_date, sim_date, 360)

        mean  = x_prev * np.exp(-mean_reversion * dt)
        var   = (sigma**2 / (2 * mean_reversion)) * (1 - np.exp(-2 * mean_reversion * dt))

        x_current = mean + np.sqrt(var) * Z[:, i]

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
    needed_dates = pd.DatetimeIndex(list(pricing_grid) + [valuation_date]).unique()
    missing_dates = [d for d in needed_dates if d not in discount_factors.index]
    if missing_dates:
        interp_vals = [
            get_discount_factor_by_zero_rates_linear_interp(
                reference_date,
                d,    
                discount_factors.index,
                discount_factors.values
                ) for d in missing_dates
            ]


        interpolated_dfs = pd.Series(index=missing_dates, data=interp_vals)

        discount_factors = (
            pd.concat([discount_factors, interpolated_dfs])
            .sort_index()
        )

    # Volatility integrand
    sigma_func = lambda u, T: (sigma / mean_reversion) * (1 - np.exp(-mean_reversion * (T - u)))    #MISTAKE? there was only one argument

    # Compute A(t,T) and C(t,T)
    for T in pricing_grid:

        if T <= valuation_date:
            A[T] = 1
            C[T] = 0
            continue

        tau_t_T = year_frac_act_x(valuation_date, T, 360)
        tau_0_T = year_frac_act_x(reference_date, T, 360)
        tau_0_t = year_frac_act_x(reference_date, valuation_date, 360)

        # integral of sigma(u,T)^2 - sigma(u,t)^2 from 0 to tau_0_t
        integral, _ = quad(lambda u: sigma_func(u, tau_0_T)**2 - sigma_func(u, tau_0_t)**2, 0, tau_0_t)
        A[T] = (discount_factors.loc[T] / discount_factors.loc[valuation_date]) * np.exp(-0.5 * integral)
        C[T] = (1 - np.exp(-mean_reversion * tau_t_T)) / mean_reversion

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
    
    # Compute target VM* (before any constraint)
    # VM* = max(MtM - THR_C, 0) if MtM >= 0
    #       min(MtM + THR_B, 0) if MtM < 0
    VM_star = np.where(mtm >= 0,
                       np.maximum(mtm - THR_C, 0),
                       np.minimum(mtm + THR_B, 0))

    # Theoretical margin call
    delta_VM = VM_star - VM

    # Effective margin call (compare with MTA)
    # delta_VM > 0: C must post -> compare with MTA_C
    # delta_VM < 0: B must post -> compare with MTA_B
    margin_call = np.where(delta_VM >= 0,
                           np.where(delta_VM >= MTA_C, delta_VM, 0),
                           np.where(np.abs(delta_VM) >= MTA_B, delta_VM, 0))

    # Update VM before caps
    VM_no_cap = VM + margin_call

    # Apply collateral caps
    # VM > 0: C has posted collateral -> cap at Cap_C
    # VM < 0: B has posted collateral -> cap at Cap_B (in absolute value)
    VM_updated = np.where(VM_no_cap >= 0,
                          np.minimum(VM_no_cap, Cap_C),
                          np.maximum(VM_no_cap, -Cap_B))

    return VM_updated



def compute_epe_base(discount_factors_input, sigma_input, simulation_grid, x_paths, fixed_leg_payment_dates, notional, mean_reversion, swap_rate_fixed):
    """
    Recompute the EPE profile for the base (10y payer) swap only,
    given a (possibly bumped) discount-factor curve and/or sigma.
    x_paths is kept fixed (same random seed) so that only the
    sensitivity to the bumped parameter is captured.
    """
    epe = np.zeros(len(simulation_grid))

    # Recompute fixed-leg cash flows on the (possibly bumped) curve
    swap_rate_bumped = swap_rate_bumped = swap_rate_fixed  # fixed contractually, does not change with the curve
    cf_bumped = pd.Series(
        data=[
            swap_rate_bumped * year_frac_30e_360(
                fixed_leg_payment_dates[k - 1], fixed_leg_payment_dates[k]
            )
            for k in range(1, len(fixed_leg_payment_dates))
        ],
        index=fixed_leg_payment_dates[1:],
    )
    cf_bumped.iloc[-1] += 1.0  # notional repayment at maturity
    cf_bumped_ext = cf_bumped.reindex(simulation_grid, fill_value=0.0)

    for i, sim_date in enumerate(simulation_grid):
        x_t = x_paths[:, i]

        # Affine components with (possibly bumped) curve and/or sigma
        A, C = affine_trick(
            sim_date, simulation_grid, mean_reversion, sigma_input, discount_factors_input
        )

        df_hw = A.values[np.newaxis, :] * np.exp(-np.outer(x_t, C.values))

        # MtM of the 10y payer swap (floating leg = par)
        dcf = df_hw[:, i:] @ cf_bumped_ext.values[i:]
        mtf = notional * (1.0 - dcf)

        epe[i] = np.mean(np.maximum(mtf, 0.0))

    return epe


def compute_cva(epe_array, lambda_input, discount_factors_input, simulation_grid, settlement_date, lgd):
    """
    compute CVA given an EPE array and a hazard rate
    Discrete CVA integral:
      CVA = Σ_i  LGD · EPE_i · B(0,t_i) · ΔPD_i
    """
    cva = 0.0
    for i, sim_date in enumerate(simulation_grid):
        prev_date = settlement_date if i == 0 else simulation_grid[i - 1]

        t_prev = year_frac_act_x(settlement_date, prev_date, 360)
        t_curr = year_frac_act_x(settlement_date, sim_date, 360)

        marginal_pd = np.exp(-lambda_input * t_prev) - np.exp(-lambda_input * t_curr)

        df_t = get_discount_factor_by_zero_rates_linear_interp(
            settlement_date,
            sim_date,
            discount_factors_input.index,
            discount_factors_input.values,
        )

        cva += lgd * epe_array[i] * df_t * marginal_pd

    return cva