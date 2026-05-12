"""
Mathematical Engineering - Financial Engineering, FY 2025-2026
Risk Management - Exercise 0: Discount Factors Bootstrap
"""

from concurrent import futures

import numpy as np
import pandas as pd
import datetime as dt
import QuantLib as ql
from utilities.date_functions import (
    business_date_offset,
    year_frac_act_x,
    year_frac_30e_360
)
from typing import Iterable, Union, List, Union, Tuple

def from_discount_factors_to_zero_rates(
    dates: Union[List[float], pd.DatetimeIndex],
    discount_factors: Iterable[float],
    convention: str = 'act/360'
) -> List[float]: 
    """
    Compute the zero rates from the discount factors.

    Parameters:
        dates (Union[List[float], pd.DatetimeIndex]): List of year fractions or dates.
        discount_factors (Iterable[float]): List of discount factors.
        convention (str): Day count convention, e.g., 'act/360' or '30/360'. Default is 'act/360'.

    Returns:
        List[float]: List of zero rates.
    """

    effDates, effDf = dates, discount_factors
    if isinstance(effDates, pd.DatetimeIndex):   # if the input are dates, convert to year fractions using the custom functions
        reference_date = effDates[0]
        
        effDates = [year_frac_act_x(reference_date, d, 360) if convention == 'act/360' else year_frac_30e_360(reference_date, d) for d in effDates[1:]]
        effDf = discount_factors[1:]
    else:  # already year fractions
        effDates = effDates[1:]
        effDf = discount_factors[1:]

    # Compute zero rates
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
    interpolation.

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
    
    
    # compute relevant yearfractions for available set of dates
    T = np.array([year_frac_act_x(reference_date, d, 365) for d in dates])
    T_target = year_frac_act_x(reference_date, interp_date, 365)    
    # convert discounts into zero rates
    zero_rates = [-np.log(df) / yf for df, yf in zip(discount_factors, T)] 
    # apply the interpolation on the target day
    interp_rate = np.interp(T_target, T, zero_rates)    
    # convert zero rate into discount
    discount = np.exp(-interp_rate * T_target)
    
    return discount 


def bootstrap(
    reference_date: dt.datetime,
    depo: pd.DataFrame,
    futures: pd.DataFrame,
    swaps: pd.DataFrame,
    shock: float = 0.0,
) -> pd.Series:
    """
    Bootstrap the discount factors from the given bid/ask market data. Deposit rates are used until
    the first future settlement date (included), futures rates are used until the 2y-swap settlement.

    Parameters:
        reference_date (dt.datetime): Reference date.
        depo (pd.DataFrame): Deposit rates.
        futures (pd.DataFrame): Futures rates.
        swaps (pd.DataFrame): Swaps rates.
        shock (Union[float, pd.Series]): Parallel shift to apply to the market rates, default to
            zero.

    Returns:
        pd.Series: Discount factors.
    """

    # initialize the list of dates and discounts
    termDates, discounts = [reference_date], [1.0]

    #### DEPOS
    
    # select the correct depos and their rates
    depoDates = depo.iloc[0:3].index.to_list()    
    depoRates = depo.loc[depoDates].mean(axis=1).values 

    # needed for the bumped bootstrap: if shock is a float, shift all the mkt data by that number, otherwise for each pillar its value
    depoRates = depoRates +  ( shock if isinstance(shock, float) else shock[depoDates].values )

    # compute B(t0, ti) = 1 / (1 + delta * L) for each deposit pillar
    termDates += depoDates
    for i in range(len(depoDates)):
        t_i  = depoDates[i]
        L_i  = depoRates[i]
        tau_i = year_frac_act_x(reference_date, t_i, 360)  # Act/360 year fraction
        B_i  = 1.0 / (1.0 + L_i * tau_i)
        discounts.append(float(B_i))


   #### FUTURES

    # select the first 7 contracts (most liquid)
    futures_of_interest = futures.iloc[0:7, :]

    # compute mid price and extract forward rate L(t0; ti-1, ti) = (100 - price) / 100
    price = (futures_of_interest["ASK"] + futures_of_interest["BID"]) / 2
    L = (100 - price) / 100

    # convert forward rate to forward discount B(t0; ti-1, ti) = 1 / (1 + L * 0.25) 
    # is a 3-month forward rate, so we use 0.25 as year fraction
    B = 1 / (1 + L * 0.25)

    # chain forward discounts to obtain spot discounts B(t0, ti) = B(t0, ti-1) * B(t0; ti-1, ti)
    for i in range(len(B)):
        discounts.append(B.iloc[i] * discounts[-1])

    # append end-of-period dates (settlement + 3 months + 2 days)
    termDates += [(d + pd.DateOffset(months=3) + pd.Timedelta(days=2)) for d in futures_of_interest.index]



    #### SWAPS

    # initialize BPV accumulator using the 1y swap date (skipped in bootstrap, covered by futures)
    swapDate = swaps.index[0]
    
    swapYearFrac, swapDisc, df = list(), list(), 0.0

    swapYearFrac.append(year_frac_30e_360(reference_date, swapDate))
    
    swap_old = swapDate
    swapDisc=[df]

    # we don't consider the first swap since we already have futures (more liquid) to cover this date
    swapRates = swaps.iloc[1:, :].mean(axis=1).values + (
        shock if isinstance(shock, float) else shock[swaps.index].values
    )
    
    # initialize BPV_1 = delta(t0, t1) * B(t0, t1) by interpolation on the existing curve
    BPV_1 = swapYearFrac[0] * get_discount_factor_by_zero_rates_linear_interp(
        reference_date, swapDate, termDates, discounts
    )
    

    # we don't consider the first swap since we already have futures (more liquid) to cover this date
    for idx, swapDate in enumerate(swaps.iloc[1:, :].index):
            rate, yf = swapRates[idx], year_frac_30e_360(swap_old, swapDate)
            swapYearFrac.append(yf)
            df = (1 - rate/100*BPV_1) / (1 + rate/100*yf)
            termDates.append(swapDate)
            discounts.append(df)
            swap_old = swapDate
            BPV_1 += df*yf

   

    discount_factors = pd.Series(index=termDates, data=discounts)
    #print('discount_factors', discount_factors)
    zero = from_discount_factors_to_zero_rates(discount_factors.index, discount_factors.values)
    zero_rates = pd.Series(index=termDates[1:], data=zero)
    return discount_factors, zero_rates
