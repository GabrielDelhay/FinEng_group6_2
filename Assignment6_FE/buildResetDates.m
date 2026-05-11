function resetDates = buildResetDates(t0, N)
% Builds the vector of N+1 quarterly reset dates starting from t0,
% adjusted to the next business day (modified following convention).

% Inputs: t0 (settlement date)
%          N (number of periods)

    resetDates    = zeros(1, N+1);
    resetDates(1) = t0;
    for i = 1:N
        resetDates(i+1) = busdate(addtodate(t0, 3*i, 'month'), 'modifiedfollow');
    end
end
