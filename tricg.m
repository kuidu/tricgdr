function [x, y, exitflag, resvec] = tricg(A, b, c, M, N, tol, maxit, opt)
% An implementation for `TRICG` in the paper:
%   "tricg and trimr: two iterative methods for symmetric qusi-definite systems"
%
% This method solve the following system:
%       [M     A] [x]   [b]
%       [       ] [ ] = [ ]
%       [A'   -N] [y]   [c]
%
% Syntaxes
% --------------
% [x, y] = tricg(A, b, c)
% [x, y] = tricg(A, b, c, M, N)
% [x, y] = tricg(A, b, c, M, N, tol)
% [x, y] = tricg(A, b, c, M, N, tol, maxit)
% [x, y] = tricg(A, b, c, M, N, tol, maxit, mreorth)
%
%
% [x, y, exitflag] = tricg(__)
% [x, y, exitflag, resvec] = tricg(__)
%
%
% Parameters
% ------------
% A                 m-by-n metrix or a function. If A is a funtion,
%                   A(x, 1) = A * x; A(x, 2) = A' * x;
%
% M, N              M and N are m-by-m and n-by-n symmetric positive matrices
%                   or funtions; If M is funtions, M(x, 1) = M * x;
%                   M(x, 2) = M \ x; N behaves in a same way, if N is funtion.
%
% b, c              m-by-1 and n-by-1 vector, respectively
%
% tol               stop tolerance
%
% maxit             maximun iterations
%
% mreorth           number of basis vectors reorthogonalized
%
% exres             determine whether compute exact residual norms or not
%
%
% Returns
%-------------
% x, y          approximate solution for the above system
%
% exitflag      covergency flag, `1` means failed, `0` means successed
%
% resvec        vectors formed by norm of residuals

arguments
    A
    b           (:, 1)          {mustBeVector}
    c           (:, 1)          {mustBeVector}
    M                                                           =  []
    N                                                           =  []
    tol         (1, 1)  double  {mustBeReal, mustBePositive}    =  1e-6
    maxit       (1, 1)  double  {mustBeInteger}                 =  100
    opt.augU                                                    =  []
    opt.augMU                                                   =  []
    opt.augV                                                    =  []
    opt.augNV                                                   =  []
    opt.mreorth (1, 1)  double  {mustBeInteger}                 =  0
    opt.exres   (1, 1)  logical                                 =  false
    opt.restart (1, 1)  double  {mustBeInteger, mustBePositive} =  1
    opt.show_info       logical                                 =  false
end

mreorth = opt.mreorth;
exres = opt.exres;
outer = opt.restart;
inner = maxit;
show_info = opt.show_info;

augU = opt.augU;
augMU = opt.augMU;
augV = opt.augV;
augNV = opt.augNV;

if isempty(M)
    M = @(x, t) x;
end

if isempty(N)
    N = @(x, t) x;
end

MisFunc = isa(M, 'function_handle');
NisFunc = isa(N, 'function_handle');
AisFunc = isa(A, 'function_handle');

x = zeros(size(b));
y = zeros(size(c));

m = length(b);
n = length(c);

if mreorth > 0
    V = zeros(m, mreorth);
    U = zeros(n, mreorth);
    MV = zeros(m, mreorth);
    NU = zeros(n, mreorth);
end

resvec = zeros(outer*inner+1, 1);
idx_res = 0;
% if exres
nr = exact_norm(A, M, N, b, c, x, y);
% else
    % nr = hypot(beta1, gamma1);
% end

idx_res = idx_res + 1;
resvec(idx_res) = nr;


exitflag = 1;

for outerit = 1:outer
    % initialization
    gx_1 = zeros(m, 1);
    gx0 = gx_1;
    gy_1 = zeros(n, 1);
    gy0 = gy_1;

    Mv0 = gx0;
    Nu0 = gy0;

    if MisFunc; Mx = M(x, 1); else; Mx = M *x; end
    if NisFunc; Ny = N(y, 1); else; Ny = N * y; end
    if AisFunc
        Ay = A(y, 1);
        ATx = A(x, 2);
    else
        Ay = A * y;
        ATx = A' * x;
    end

    Mv1 = b - (Mx + Ay);
    Nu1 = c - (ATx - Ny);

    if MisFunc; v1 = M(Mv1, 2); else; v1 = M \ Mv1; end
    if NisFunc; u1 = N(Nu1, 2); else; u1 = N \ Nu1; end

    beta1 = sqrt(v1' * Mv1);
    gamma1 = sqrt(u1' * Nu1);

    Mv1 = Mv1 / beta1;
    Nu1 = Nu1 / gamma1;
    v1 = v1 / beta1;
    u1 = u1 / gamma1;

    d_1 = 0;
    d0 = 0;
    sig1 = 0;
    eta1 = 0;
    lam1 = 0;
    delta0 = 0;
    pi_1 = 0;
    pi0 = 0;


    last_col = 0;

    for k = 1:inner

        if mreorth > 0
            idx = mod(k-1, mreorth) + 1;
            MV(:, idx) = Mv1;
            NU(:, idx) = Nu1;
            V(:, idx) = v1;
            U(:, idx) = u1;

            last_col = min(last_col + 1, mreorth);
        end

        if AisFunc
            Au = A(u1, 1);
            ATv = A(v1, 2);
        else
            Au = A * u1;
            ATv = A' * v1;
        end

        q = Au - gamma1 * Mv0;
        p = ATv - beta1 * Nu0;

        alpha = v1' * q;

        Mv = q - alpha * Mv1;
        Nu = p - alpha' * Nu1;

        if mreorth > 0
            Mv = reorth(Mv, V(:, 1:last_col), MV(:, 1:last_col));
            Nu = reorth(Nu, U(:, 1:last_col), NU(:, 1:last_col));
        end

        if ~isempty(augU) && ~isempty(augMU)
            Mv = Mv - augMU * (augU'*Mv);
        end

        if ~isempty(augV) && ~isempty(augNV)
            Nu = Nu - augNV * (augV'*Nu);
        end
            
        if MisFunc; v = M(Mv, 2); else; v = M \ Mv; end
        if NisFunc; u = N(Nu, 2); else; u = N \ Nu; end

        beta2 = sqrt(v' * Mv);
        gamma2 = sqrt(u' * Nu);

        if k == 1
            d1 = 1;
            delta1 = alpha' / d1;
            d2 = -1 - abs(delta1)^2 * d1;
        else
            sig1 = beta1 / d0;
            eta1 = gamma1 / d_1;
            lam1 = -(eta1 * delta0' * d_1) / d0;
            d1 = 1 - abs(sig1)^2 * d0;
            delta1 = (alpha' - lam1 * sig1' * d0) / d1;
            d2 = -1 - abs(eta1)^2 * d_1 - abs(lam1)^2 * d0 - abs(delta1)^2 * d1;
        end

        if k == 1
            pi1 = beta1 / d1;
            pi2 = (gamma1 - delta1 * beta1) / d2;
        else
            pi1 = -(sig1 * pi0 * d0) / d1;
            pi2 = -(delta1 * pi1 * d1 + lam1 * pi0 * d0 + eta1 * pi_1 * d_1) / d2;
        end

        % update approximate solution
        gx1 = v1 - sig1' * gx0;
        gx2 = -delta1' * gx1 - lam1' * gx0 - eta1' * gx_1;

        gy1 = -sig1' * gy0;
        gy2 = u1 - delta1' * gy1 - lam1' * gy0 - eta1' * gy_1;

        x = x + pi1 * gx1 + pi2 * gx2;
        y = y + pi1 * gy1 + pi2 * gy2;

        % compute norm of residual exactly
        if exres
            nr = exact_norm(A, M, N, b, c, x, y);
        else
            nr = hypot(gamma2 * (pi1 - delta1' * pi2), beta2 * pi2);
        end
        idx_res = idx_res + 1;
        resvec(idx_res) = nr;

        if nr <= tol
            if ~exres
                nr = exact_norm(A, M, N, b, c, x, y);
                resvec(idx_res) = nr;
            end

            if nr <= tol
                exitflag = 0;
                resvec = resvec(1:idx_res);
                return;
            end
        end

        % update for next iteration
        Mv0 = Mv1;
        Nu0 = Nu1;
        Mv1 = Mv / beta2;
        Nu1 = Nu / gamma2;
        v1 = v / beta2;
        u1 = u / gamma2;
        gamma1 = gamma2;
        beta1 = beta2;

        d_1 = d1; d0 = d2;
        pi_1 = pi1; pi0 = pi2;
        delta0 = delta1;

        gx_1 = gx1; gx0 = gx2;
        gy_1 = gy1; gy0 = gy2;
    end
end

if show_info
    if exitflag == 1
        fprintf(['Reach maximum iteration: %d. No covergence to required ' ...
            'tolerance %e. Absolute error is %.4e.\n'], outerit*inner, tol, resvec(idx_res));
    end
end

% resvec = resvec(1:idx_res);

end