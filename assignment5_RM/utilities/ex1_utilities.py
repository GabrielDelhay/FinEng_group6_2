"""
Mathematical Engineering - Financial Engineering, FY 2025-2026
Risk Management - Exercise 1: Hedging a Swaption Portfolio
"""

from __future__ import annotations

from enum import Enum
from typing import List, Optional, Tuple, Union

import datetime as dt

import numpy as np
import pandas as pd
from scipy.stats import norm

from utilities.date_functions import (
    date_series,
    schedule_year_fraction,
    year_frac_30e_360,
    year_frac_act_x,
)
from utilities.ex0_utilities import get_discount_factor_by_zero_rates_linear_interp

__all__ = [
    "SwapType",
    "basis_point_value",
    "swap_par_rate",
    "swap_mtm",
    "irs_proxy_duration",
    "swaption_price_calculator",
]

# Discount factors
def _df(
    discount_factors: pd.Series,
    date: Union[dt.date, pd.Timestamp],
) -> float:
    """Thin wrapper around the linear-interpolation discount-factor lookup."""
    return get_discount_factor_by_zero_rates_linear_interp(
        discount_factors.index[0],
        date,
        discount_factors.index,
        discount_factors.values,
    )


# Enumerations
class SwapType(Enum):
    RECEIVER = "receiver"
    PAYER = "payer"



# Core building blocks

def basis_point_value(
    fixed_leg_schedule: List[Union[dt.datetime, dt.date]],
    discount_factors: pd.Series,
    settlement_date: Optional[Union[dt.datetime, dt.date]] = None,
) -> float:
    """
    Annuity / basis-point value (BPV) of a swap fixed leg.

    The BPV is defined as:

        A(t, T) = sum_{i=1}^{N} tau_i * P(t, T_i)

    where tau_i = year_frac_30E_360(T_{i-1}, T_i) and P(t, T_i) is the
    discount factor seen from the value date t.

    When *settlement_date* is provided (forward-starting swap), all discount
    factors are divided by P(t, T_settlement) to express the BPV as seen from
    the forward settlement date.

    Parameters
    ----------
    fixed_leg_schedule:
        Coupon payment dates T_1, …, T_N (i.e. excluding the start date T_0).
    discount_factors:
        Zero-curve discount factors indexed by date.
    settlement_date:
        Forward settlement date T_0.  If None the value date
        ``discount_factors.index[0]`` is used.

    Returns
    -------
    float
        Annuity factor (BPV).
    """
    if len(fixed_leg_schedule) == 0:
        raise ValueError("fixed_leg_schedule must contain at least one date.")

    # Reference date for the day-count start of the first coupon period
    if settlement_date is None:
        t0 = discount_factors.index[0]
        df_t0 = 1.0  # no forward adjustment needed
    else:
        t0 = settlement_date
        df_t0 = _df(discount_factors, t0)
        if df_t0 <= 0:
            raise ValueError(f"Discount factor at settlement_date {t0} is non-positive.")

    # Build the full date sequence: [T_0, T_1, …, T_N]
    dates = [t0] + list(fixed_leg_schedule)

    bpv = sum(
        year_frac_30e_360(dates[i], dates[i + 1]) * _df(discount_factors, dates[i + 1])
        for i in range(len(dates) - 1)
    ) / df_t0

    return bpv


def swap_par_rate(
    fixed_leg_schedule: List[Union[dt.datetime, dt.date]],
    discount_factors: pd.Series,
    fwd_start_date: Optional[Union[dt.datetime, dt.date]] = None,
) -> float:
    """
    Par swap rate (or forward swap rate when *fwd_start_date* is provided).

        S = (P(t, T_0) - P(t, T_N)) / A(t, T)

    Parameters
    ----------
    fixed_leg_schedule:
        Coupon payment dates T_1, …, T_N.
    discount_factors:
        Zero-curve discount factors indexed by date.
    fwd_start_date:
        Forward start date T_0.  If None, T_0 defaults to the value date.

    Returns
    -------
    float
        Par (or forward) swap rate.
    """
    df_tN = _df(discount_factors, fixed_leg_schedule[-1])

    if fwd_start_date is not None:
        df_t0 = _df(discount_factors, fwd_start_date)
        float_leg = 1.0 - df_tN / df_t0
        bpv = basis_point_value(fixed_leg_schedule, discount_factors, fwd_start_date)
    else:
        float_leg = 1.0 - df_tN
        bpv = basis_point_value(fixed_leg_schedule, discount_factors)

    if bpv == 0:
        raise ZeroDivisionError("BPV is zero — cannot compute par rate.")

    return float_leg / bpv


def swap_mtm(
    swap_rate: float,
    payments_schedule: List[Union[dt.datetime, dt.date]],
    discount_factors: pd.Series,
    swap_type: SwapType = SwapType.PAYER,
) -> float:
    """
    Mark-to-market of a vanilla interest rate swap (single-curve framework).

    MtM = sign * (float_leg - fixed_leg)
        = sign * ((1 - P(t, T_N)) - K * A(t, T))

    where sign = +1 for payer, -1 for receiver.

    Parameters
    ----------
    swap_rate:
        Fixed coupon rate K.
    payments_schedule:
        Fixed leg payment dates T_1, …, T_N.
    discount_factors:
        Zero-curve discount factors indexed by date.
    swap_type:
        ``SwapType.PAYER`` or ``SwapType.RECEIVER``.

    Returns
    -------
    float
        Swap MtM (per unit notional).
    """
    bpv = basis_point_value(payments_schedule, discount_factors)
    p_tN = _df(discount_factors, payments_schedule[-1])
    float_leg = 1.0 - p_tN
    fixed_leg = swap_rate * bpv

    if swap_type == SwapType.PAYER:
        return float_leg - fixed_leg
    elif swap_type == SwapType.RECEIVER:
        return fixed_leg - float_leg
    else:
        raise ValueError(f"Unknown swap type: {swap_type!r}")


def irs_proxy_duration(
    ref_date: dt.date,
    swap_rate: float,
    fixed_leg_payment_dates: List[dt.date],
    discount_factors: pd.Series,
) -> float:
    """
    Rate sensitivity of an IRS, proxied by the Macaulay duration of an
    equivalent fixed-coupon bond (30E/360 day count).

    Parameters
    ----------
    ref_date:
        Valuation date.
    swap_rate:
        Fixed coupon rate of the swap.
    fixed_leg_payment_dates:
        Coupon (and final principal) payment dates.
    discount_factors:
        Zero-curve discount factors indexed by date.

    Returns
    -------
    float
        Macaulay duration (in years).
    """
    if len(fixed_leg_payment_dates) == 0:
        raise ValueError("fixed_leg_payment_dates must not be empty.")

    weighted_time = 0.0
    price = 0.0

    for i, payment_date in enumerate(fixed_leg_payment_dates):
        t = year_frac_30e_360(ref_date, payment_date)
        df = _df(discount_factors, payment_date)
        cf = swap_rate + (1.0 if i == len(fixed_leg_payment_dates) - 1 else 0.0)
        pv = cf * df
        weighted_time += t * pv
        price += pv

    if price == 0:
        raise ZeroDivisionError("Bond price is zero — cannot compute duration.")

    return weighted_time / price


def swaption_price_calculator(
    S0: float,
    strike: float,
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    underlying_expiry: Union[dt.date, pd.Timestamp],
    sigma_black: float,
    freq: int,
    discount_factors: pd.Series,
    swaption_type: SwapType = SwapType.RECEIVER,
    compute_delta: bool = False,
) -> Union[float, Tuple[float, float]]:
    """
    Black swaption price.

    Parameters
    ----------
    S0:
        Forward swap rate (annualised).
    strike:
        Fixed rate of the underlying forward swap (swaption strike).
    ref_date:
        Valuation date.
    expiry:
        Option expiry date (= swap start date).
    underlying_expiry:
        Underlying swap maturity date.
    sigma_black:
        Lognormal (Black) implied volatility.
    freq:
        Number of fixed-leg coupon payments per year.
    discount_factors:
        Zero-curve discount factors indexed by date.
    swaption_type:
        ``SwapType.PAYER`` or ``SwapType.RECEIVER``.
    compute_delta:
        If True, also return the Black delta.

    Returns
    -------
    float or (float, float)
        Swaption price in present-value terms (and delta if requested).

    Delta is taken with respect to S0:
        Payer delta   = P(t,T_exp) * A * N(d1)
        Receiver delta= P(t,T_exp) * A * (N(d1) - 1)
    """
    if S0 <= 0:
        raise ValueError(f"Forward swap rate S0 must be positive, got {S0}.")
    if strike <= 0:
        raise ValueError(f"Strike must be positive, got {strike}.")
    if sigma_black <= 0:
        raise ValueError(f"Volatility sigma_black must be positive, got {sigma_black}.")

    ttm = year_frac_act_x(ref_date, expiry, 365)
    if ttm <= 0:
        raise ValueError("Expiry must be strictly after ref_date.")

    sqrt_ttm = np.sqrt(ttm)
    d1 = (np.log(S0 / strike) + 0.5 * sigma_black ** 2 * ttm) / (sigma_black * sqrt_ttm)
    d2 = d1 - sigma_black * sqrt_ttm

    fixed_leg_dates = date_series(expiry, underlying_expiry, freq)
    # fixed_leg_dates[0] is the swap start (= expiry); coupon dates start at [1]
    annuity = basis_point_value(fixed_leg_dates[1:], discount_factors, expiry)
    p_expiry = _df(discount_factors, expiry)

    if swaption_type == SwapType.PAYER:
        price = p_expiry * annuity * (S0 * norm.cdf(d1) - strike * norm.cdf(d2))
        delta = p_expiry * annuity * norm.cdf(d1)
    elif swaption_type == SwapType.RECEIVER:
        price = p_expiry * annuity * (strike * norm.cdf(-d2) - S0 * norm.cdf(-d1))
        delta = p_expiry * annuity * (norm.cdf(d1) - 1.0)
    else:
        raise ValueError(f"Unknown swaption type: {swaption_type!r}")

    return (price, delta) if compute_delta else price

