"""
Mathematical Engineering - Financial Engineering, FY 2025-2026
Risk Management - Exercise 0: Discount Factors Bootstrap
"""

import numpy as np
import pandas as pd
import datetime as dt
from utilities.date_functions import (
    business_date_offset,
    year_frac_act_x,
    year_frac_30e_360
)
# modification: added Tuple, removed duplicate Union and unused QuantLib/concurrent imports
from typing import Iterable, Union, List, Tuple


def from_discount_factors_to_zero_rates(
    dates: Union[List[float], pd.DatetimeIndex],
    discount_factors: Iterable[float],
    # modification: default convention changed to act/365 as required by the professor
    convention: str = 'act/365'
) -> List[float]:
    """
    Compute the zero rates from the discount factors.

    Parameters:
        dates (Union[List[float], pd.DatetimeIndex]): List of year fractions or dates.
        discount_factors (Iterable[float]): List of discount factors.
        convention (str): Day count convention, 'act/365' or 'act/360'. Default is 'act/365'.

    Returns:
        List[float]: List of zero rates.
    """
    effDates, effDf = dates, discount_factors
    if isinstance(effDates, pd.DatetimeIndex):
        reference_date = effDates[0]
        # modification: act/365 convention applied consistently (act/360 was incorrect)
        effDates = [
            year_frac_act_x(reference_date, d, 365) if convention == 'act/365'
            else year_frac_act_x(reference_date, d, 360)
            for d in effDates[1:]
        ]
        effDf = discount_factors[1:]
    else:
        effDates = effDates[1:]
        effDf = discount_factors[1:]

    zero_rates = [-np.log(df) / yf for df, yf in zip(effDf, effDates)]
    return zero_rates


def get_discount_factor_by_zero_rates_linear_interp(
    reference_date: Union[dt.datetime, pd.Timestamp],
    interp_date: Union[dt.datetime, pd.Timestamp],
    dates: Union[List[dt.datetime], pd.DatetimeIndex],
    discount_factors: Iterable[float],
) -> float:
    """
    Given a list of discount factors, return the discount factor at a given date by linear
    interpolation on zero rates (Act/365).

    Parameters:
        reference_date (Union[dt.datetime, pd.Timestamp]): Reference date.
        interp_date (Union[dt.datetime, pd.Timestamp]): Date at which the discount factor is
            interpolated.
        dates (Union[List[dt.datetime], pd.DatetimeIndex]): List of dates.
        discount_factors (Iterable[float]): List of discount factors.

    Returns:
        float: Discount factor at the interpolated date.
    """
    if len(dates) != len(discount_factors):
        raise ValueError("Dates and discount factors must have the same length.")

    dates_index = pd.DatetimeIndex(dates)
    T_target    = year_frac_act_x(reference_date, interp_date, 365)

    # modification: reuse from_discount_factors_to_zero_rates instead of duplicating
    # the discount-to-zero-rate conversion logic
    zero_rates_vals = from_discount_factors_to_zero_rates(
        dates_index, np.array(discount_factors), convention='act/365'
    )

    # modification: prepend t0 point (T=0, zero rate=0) to allow interpolation from
    # reference date — without this, np.interp cannot handle dates before the first pillar
    T_all     = np.array([year_frac_act_x(reference_date, d, 365) for d in dates_index])
    T_interp  = np.concatenate([[0.0], T_all[1:]])
    zr_interp = np.concatenate([[0.0], zero_rates_vals])

    interp_rate = np.interp(T_target, T_interp, zr_interp)
    return np.exp(-interp_rate * T_target)


def bootstrap(
    reference_date: dt.datetime,
    depo: pd.DataFrame,
    futures: pd.DataFrame,
    swaps: pd.DataFrame,
    shock: float = 0.0,
    # modification: return type corrected from pd.Series to Tuple[pd.Series, pd.Series]
) -> Tuple[pd.Series, pd.Series]:
    """
    Bootstrap the discount factors from the given bid/ask market data.
    Deposits cover up to and including the first future settlement date,
    first 7 futures cover the mid range, swaps (2y onward) cover the long end.

    Parameters:
        reference_date (dt.datetime): Reference date.
        depo (pd.DataFrame): Deposit rates (decimal), BID/ASK columns, index = maturity dates.
        futures (pd.DataFrame): Futures joined with expiry, BID/ASK/Settle/Expiry columns,
            index = IMM settlement dates.
        swaps (pd.DataFrame): Swap rates (percent), BID/ASK columns, index = maturity dates.
        shock (Union[float, pd.Series]): Parallel shift (decimal) applied to all rates.

    Returns:
        Tuple[pd.Series, pd.Series]: Discount factors and zero rates indexed by date.
    """

    termDates = [reference_date]
    discounts = [1.0]

    # -------------------------------------------------------------------------
    # DEPOSITS  B(t0, Ti) = 1 / (1 + L * tau),  tau = Act/360
    # -------------------------------------------------------------------------

    # modification: dynamic selection — deposits maturing up to and including the first
    # future settlement date, avoids hard-coded iloc[0:3]
    first_future_settle = futures["Settle"].iloc[0]
    depo_selected       = depo[depo.index <= first_future_settle]

    depoDates = depo_selected.index.to_list()
    depoRates = depo_selected.mean(axis=1).values

    # apply parallel shock to deposit rates (already in decimal)
    depoRates = depoRates + (shock if isinstance(shock, float) else shock[depoDates].values)

    # modification: vectorised computation — for loop replaced by array operations
    taus     = np.array([year_frac_act_x(reference_date, t, 360) for t in depoDates])
    depo_dfs = 1.0 / (1.0 + depoRates * taus)

    termDates += depoDates
    discounts += depo_dfs.tolist()

    # modification: if first future settlement date is not already a deposit maturity,
    # add it explicitly by interpolation — required output point, boundary between segments
    if first_future_settle not in depoDates:
        B_first_settle = get_discount_factor_by_zero_rates_linear_interp(
            reference_date, first_future_settle, termDates, discounts
        )
        termDates.append(first_future_settle)
        discounts.append(B_first_settle)

    # -------------------------------------------------------------------------
    # STIR FUTURES  B(t0, T_expiry) = B(t0, T_settle) / (1 + L * tau),  tau = Act/360
    # -------------------------------------------------------------------------

    # modification: fixed selection of the first 7 futures contracts as required
    futures_of_interest = futures.iloc[0:7]
    settle_dates        = futures_of_interest["Settle"]
    expiry_dates        = futures_of_interest["Expiry"]

    price = (futures_of_interest["ASK"] + futures_of_interest["BID"]) / 2
    L     = (100 - price) / 100

    # apply parallel shock to futures forward rates
    L = L + (shock if isinstance(shock, float) else shock[futures_of_interest.index].values)

    # modification: sequential chaining — B(t0, expiry_i) = B(t0, expiry_{i-1}) * B_fwd_i
    # B(t0, T_settle_i) obtained by interpolation on the current curve at each step
    # modification: exact Act/360 year fraction between Settle and Expiry (not 0.25)
    for i, (settle, expiry) in enumerate(zip(settle_dates, expiry_dates)):
        B_settle = get_discount_factor_by_zero_rates_linear_interp(
            reference_date, settle, termDates, discounts
        )
        tau_i = year_frac_act_x(settle, expiry, 360)
        B_fwd = 1.0 / (1.0 + L.iloc[i] * tau_i)
        # modification: only expiry dates stored — not settlement dates
        termDates.append(expiry)
        discounts.append(B_settle * B_fwd)

    # -------------------------------------------------------------------------
    # SWAPS  B(t0, Tn) = (1 - K * BPV_{n-1}) / (1 + K * tau_n),  tau = 30/360
    # -------------------------------------------------------------------------

    # modification: swap rates kept in percent here, /100 applied inside loop
    # consistent with original code — shock must be in percent units if rates are in percent
    swapRates = swaps.iloc[1:, :].mean(axis=1).values + (
        shock if isinstance(shock, float) else shock[swaps.index[1:]].values
    )

    first_swap_date = swaps.index[0]
    swap_old        = first_swap_date

    # modification: BPV initialised with the 1y swap point — covered by futures, NOT added
    # to the output curve, but used in BPV calculation as required
    yf_1 = year_frac_30e_360(reference_date, first_swap_date)
    B_1  = get_discount_factor_by_zero_rates_linear_interp(
        reference_date, first_swap_date, termDates, discounts
    )
    BPV = yf_1 * B_1

    # modification: loop from 2nd swap onward — 1y covered by futures
    # modification: removed unused variables swapYearFrac, swapDisc, df=0.0
    for swap_date, rate in zip(swaps.index[1:], swapRates):
        yf = year_frac_30e_360(swap_old, swap_date)
        df = (1.0 - rate / 100.0 * BPV) / (1.0 + rate / 100.0 * yf)
        termDates.append(swap_date)
        discounts.append(df)
        # modification: BPV updated incrementally — BPV_N = BPV_{N-1} + tau_N * B_N
        BPV     += yf * df
        swap_old = swap_date

    # -------------------------------------------------------------------------
    # OUTPUT
    # -------------------------------------------------------------------------

    # modification: explicit sort_index() — guarantees chronological order
    discount_factors = pd.Series(index=termDates, data=discounts).sort_index()

    # modification: NaN guard — raises immediately if bootstrap produced invalid values
    if discount_factors.isna().any():
        raise ValueError("NaN detected in discount factors — check input data")

    # modification: act/365 passed explicitly to zero rate conversion
    zero = from_discount_factors_to_zero_rates(
        discount_factors.index, discount_factors.values, convention='act/365'
    )
    # modification: zero rates index aligned on sorted discount_factors index
    zero_rates = pd.Series(index=discount_factors.index[1:], data=zero)

    return discount_factors, zero_rates