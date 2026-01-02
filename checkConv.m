function [conv, absres, relres, len_conv] = checkConv(res1, res2, k, tol, Smax)
% This function checks the convergence of singular values

conv = false;
% tol = Smax * tol;

I = (res1(1:k) <= tol);
len_conv = length( find(res2(I) <= tol) );

absres = max([res1(I), res2(I)]);
relres = absres / Smax;

if len_conv == k
    % U_T = U_T(:, 1:K_org); V_T = V_T(:, 1:K_org); S_T = S_T(1:K_org);
    conv = true;
end
end

