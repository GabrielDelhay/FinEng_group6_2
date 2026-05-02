function resetDates = buildResetDates(t0, N)
    resetDates    = zeros(1, N+1);
    resetDates(1) = t0;
    for i = 1:N
        resetDates(i+1) = busdate(addtodate(t0, 3*i, 'month'), 'modifiedfollow');
    end
end
