function [errAbs, errRel] = pNormError(x1, x2, p)


errRel = norm(x1 - x2,p)/norm(x2,p); % Relative error
errAbs = norm(x1 - x2,p); % Absolute error

