function NPV_A = price_floating_leg(notional, B_pay, delta_pay, spread)
%  INPUTS:
%    notional   : N 
%    B_start    : discount factor B(t0, T_start)  
%    B_pay      : discount factors B(t0, T_i) at the N payment dates
%    delta_pay  : Act/360 year fractions  delta_i = delta(T_{i-1}, T_i)
%    spread     : annualized spread 
%
%  OUTPUT:
%    NPV        : present value of the floating leg, in currency units

floater_pv = notional * (1 - B_pay(end));
BPV        = sum(delta_pay(:) .* B_pay(:));
spread_pv  = notional * spread * BPV;

NPV_A = floater_pv + spread_pv;

end