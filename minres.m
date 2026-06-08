function [x, exitflag, resvec] = minres(A, b, M, tol, maxit, opt)
% An implementation for `MINRES`. MINRES solve the linear system
%       Ax = b
%
% or least square problems
%       min || Ax - b ||_2
%
% MINRES finds the solution x which is belong to K(A,b).
%
%
% Syntaxes
% --------------
% x = minres(A, b)
% x = minres(A, b, M)
% x = minres(A, b, M, tol)
% x = minres(A, b, M, tol, maxit)
% x = minres(A, b, M, tol, maxit, x0)
%
%
% [x, exitflag] = minres(__)
% [x, exitflag, resvec] = minres(__)
%
%
% Parameters
% ------------
% A                 m by n matrix
%
% b                 m by 1 right hand vector
%
% tol               stop tolerance
%
% maxit             number of maximum iterations
%
% Name, Value       Name-Value pairs determine other options. Name should be in
%                   {'x0'}. 'x0' means initial gauss; Corresponding default
%                   values are {[]}.
%
%
% Returns
%-------------
% x                 approximate solution for the above system
%
% exitflag          convergence flag, `1` means failed, `0` means successes
%
% resvec            vectors formed by norm of residuals

arguments
    A                                {mustBeNonempty}
    b                                {mustBeNonempty}
    M                                                                = []
    tol              (1, 1)  double  {mustBeNonnegative}             = 1e-6
    maxit            (1, 1)  double  {mustBeInteger, mustBePositive} = 100
    opt.exres        (1, 1)  logical                                 = false
    opt.x0           (:, 1)  double  {mustBeVector}                  = []
end

exres = opt.exres;

if isempty(M)
    Mfunc = @(x, t) x;
elseif isa(M, 'function_handle')
    Mfunc = M;
else
    funcs = {@(x) M*x, @(x) M\x};
    Mfunc = @(x, t) funcs{mod(t-1,2)+1}(x);
end

if isa(A, 'function_handle')
    Afunc = A;
else
    Afunc = @(x) A * x;
end

m = length(b);
x0 = opt.x0;
if isempty(x0)
    x = zeros(m, 1);
end

r = b - Afunc(x);

v1 = Mfunc(r, 2);
beta1 = sqrt(r' * v1);
v1 = v1 / beta1;
Mv1 = r / beta1;
Mv0 = zeros(m, 1);

p0 = zeros(m, 1); 
p_1 = p0;

eta_1 = 0; 
tautilde = beta1;

exitflag = 1;
resvec = zeros(maxit+1, 1);
resvec(1) = beta1;
for k = 1:maxit
    q = A * v1 - beta1 * Mv0;
    alpha = v1' * q;

    Mv = q - alpha * Mv1;
    v = Mfunc(Mv, 2);
    beta2 = sqrt(v' * Mv);

    if k == 1
        alphatilde = alpha; 
        beta2tilde = beta2;
        lam0 = 0; eta0 = 0;
    else
        alphatilde = -s * beta1tilde + c * alpha;
        lam0 = c * beta1tilde + s * alpha;

        eta0 = s * beta2; beta2tilde = c * beta2;
    end

    [G, y] = planerot([alphatilde; beta2]);
    delta = y(1); c = G(1, 1); s = G(1, 2);

    if k == 1
        p1 = v1 / delta;
    elseif k == 2
        p1 = (v1 - lam0 * p0) / delta;
    else
        p1 = (v1 - lam0 * p0 - eta_1 * p_1) / delta;
    end
    
    tau = c * tautilde; 
    tautilde = -s * tautilde;

    dx = tau * p1;
    x = x + dx;
    
    if exres
        r = b - A * x;
        nr = sqrt(r' * Mfunc(r, 2));
    else
        nr = abs(tautilde);
    end

    resvec(k+1) = nr;
    if nr <= tol
        if ~exres
            r = b - A * x;
            nr = sqrt(r' * Mfunc(r, 2));
        end

        if nr <= tol
            exitflag = 0;
            break;
        end
    end

    % update for next iteration
    v1 = v / beta2;
    Mv0 = Mv1; 
    Mv1 = Mv / beta2;

    beta1 = beta2;
    eta_1 = eta0;
    beta1tilde = beta2tilde;

    p_1 = p0; 
    p0 = p1;
end

resvec = resvec(1:k+1);

end