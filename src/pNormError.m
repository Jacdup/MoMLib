function [err] = pNormError(x1, x2, p)


err = norm(x1 - x2,p)/norm(x2,p); % Relative error

