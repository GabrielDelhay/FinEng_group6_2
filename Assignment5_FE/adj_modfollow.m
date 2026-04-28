function adj = adj_modfollow(raw_dates)
% ADJ_MODFOLLOW  Modified Following Business Day Convention (vectorized).
%
%   adj = adj_modfollow(raw_dates) takes an array of MATLAB datenums and
%   returns the adjusted dates under the Modified Following rule:
%     - if the raw date is a business day, keep it;
%     - otherwise, move to the NEXT business day;
%     - if that next business day rolls into a new month, fall back to
%       the PREVIOUS business day instead.
%
%   The output has the same size and shape as the input.
%
%   INPUT:
%     raw_dates : array of datenums (any shape)
%
%   OUTPUT:
%     adj       : array of adjusted datenums (same shape as input)

    inputSize = size(raw_dates);
    raw       = raw_dates(:);              % flatten to column

    adj  = raw;                            % default: leave untouched
    mask = ~isbusday(raw);                 % positions that need adjustment

    if any(mask)
        next_bd = busdate(raw(mask),  1);  % next business day
        prev_bd = busdate(raw(mask), -1);  % previous business day

        % If the next BD falls in a different month -> use the previous BD
        rollover           = month(next_bd) ~= month(raw(mask));
        chosen             = next_bd;
        chosen(rollover)   = prev_bd(rollover);

        adj(mask) = chosen;
    end

    adj = reshape(adj, inputSize);         % preserve original shape
end