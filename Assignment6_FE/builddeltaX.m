function [deltaX_top, deltaX_inner, deltaX_bot] = builddeltaX(dx, lMax)
% Inner rows (2 : end-1): pattern A — branches at (l+1, l, l-1)
deltaX_inner = repmat([+dx, 0, -dx], 2*lMax-1, 1);   % (2*lMax-1) x 3

% Top row (1): pattern C — branches at (l, l-1, l-2)
deltaX_top = [0, -dx, -2*dx];                        % 1 x 3

% Bottom row (end): pattern B — branches at (l+2, l+1, l)
deltaX_bot = [+2*dx, +dx, 0];                        % 1 x 3

end