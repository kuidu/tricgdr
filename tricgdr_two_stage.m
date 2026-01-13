function [x, y, exitflag, resvec, varargout] = tricgdr_two_stage(A, b, c, M, N, k, p, tol, maxcycle, opt)
% TRICGDR will find approximate solution of following equation with
% implicity deflated restart:
%
%       [M    A] [x]   [b]
%       [      ] [ ] = [ ]
%       [A'  -N] [y]   [c]
%
%
% Syntaxes
% --------------
% [x, y] = tricgdr_two_stage(A, b, c)
% [x, y] = tricgdr_two_stage(A, b, c, M, N)
% [x, y] = tricgdr_two_stage(A, b, c, M, N, k, p)
% [x, y] = tricgdr_two_stage(A, b, c, M, N, k, p, tol)
% [x, y] = tricgdr_two_stage(A, b, c, M, N, k, p, tol, maxcycle)
% [x, y] = tricgdr_two_stage(__, Name, Value)
%
%
% [x, y, exitflag] = tricgdr_two_stage(__)
% [x, y, exitflag, resvec] = tricgdr_two_stage(__) 
% [x, y, exitflag, resvec, U, T, V, MU, NV] = tricgdr_two_stage(__)
%
% Parameters
% --------------
% A             a m x n matrix
%
% b, c          right-hand vectors
%
% M, N          SPD matrices
%
% k             number of approximate augmenting singular vectors. k is a
%               vector of length 2: [k1, k2], k1 is the number of desired
%               largest singular vectors, and k2 is the number of desired
%               smallest singular vectors
%
% p             number of dimension of subspace
%
% tol           stop tolerance
% 
% maxcycle      number of maximum iterations
%
%
% Name, Value   name-Value pairs determine other options.
%                   'stol':         tolerance for approximate elliptic singular
%                                   triplets
%                   
%                   'maxit':        number of maximum iterations on the
%                                   non-restarting stage
% 
%                   'mreorth':      number of converged singular vectors used
%                                   for reorthogonalization
% 
%                   'exres':        determine whether compute the exact
%                                   residual norms
% 
%                   'reorth':       the type of reorthogonalization. `1`
%                                   means only reorthogonalize one side.
%                                   `2` means reorthogonalize both two sides.
%                   
%                   'show_info':    determine whether print the information
%                                   of the results
%
%
% Returns
% --------------
% x, y              approximate solution for the above system
%
% exitflag          covergency flag, `1` means failed, `0` means successed
%
% resvec            vectors formed by norm of residuals 
% 
% U, T, V, MU, NV   approximate elliptic singular triplets

arguments
    A
    b               (:, 1)          {mustBeVector}
    c               (:, 1)          {mustBeVector}
    M                                                                       =  []
    N                                                                       =  []
    k                       double  {mustBeVector}                          =  [10, 0]
    p               (1, 1)  double  {mustBeInteger}                         =  40
    tol             (1, 1)  double  {mustBeReal, mustBeNonnegative}         =  1e-6
    maxcycle        (1, 1)  double  {mustBeInteger}                         =  100
    opt.stol        (1, 1)  double  {mustBeReal, mustBeNonnegative}         =  tol
    opt.maxit       (1, 1)  double  {mustBeInteger}                         =  100
    opt.mreorth     (1, 1)  double  {mustBeInteger}                         =  sum(k)
    opt.exres       (1, 1)  logical                                         =  false
    opt.reorth      (1, 1)  double  {mustBeMember(opt.reorth, [0, 2])}      =  2
    opt.show_info   (1, 1)  logical                                         =  false
end

outer = maxcycle;
stol = opt.stol;
maxit = opt.maxit;
mreorth = opt.mreorth;
reorth_type = opt.reorth;
exres = opt.exres;
show_info = opt.show_info;

if mreorth > p
    mreorth = sum(k);
end

if isempty(M)
    M = @(x, t) x;
end

if isempty(N)
    N = @(x, t) x;
end

MisFunc = isa(M, 'function_handle');
NisFunc = isa(N, 'function_handle');
AisFunc = isa(A, 'function_handle');

k1 = k(1);
k2 = k(2);

m = length(b);
n = length(c);

U = zeros(m, p);
V = zeros(n, p);
MU = zeros(m, p);
NV = zeros(n, p);
T = zeros(p, p);
x = zeros(m, 1);
y = zeros(n, 1);

Mu1 = b;
Nv1 = c;

if MisFunc; u1 = M(Mu1, 2); else; u1 = M \ Mu1; end
if NisFunc; v1 = N(Nv1, 2); else; v1 = N \ Nv1; end

beta1 = sqrt(u1' * Mu1);
gamma1 = sqrt(v1' * Nv1);

u1 = u1 / beta1;
v1 = v1 / gamma1;
Mu1 = Mu1 / beta1;
Nv1 = Nv1 / gamma1;
U(:, 1) = u1;
V(:, 1) = v1;
MU(:, 1) = Mu1;
NV(:, 1) = Nv1;

exitflag = 1;
resvec = zeros(outer*p+1, 1);

if exres
    nr = exact_norm(A, M, N, b, c, x, y);
else
    nr = hypot(beta1, gamma1);
end

resvec(1) = nr;

k_aug = k1 + k2;
k = 0;
idx_res = 1;

if k_aug >= p
    k_aug = min(max(p-10, 10), floor(p/2));
    warning(['dimension of augment space should be less than dimension of subspace. ' ...
        'set k = %d\n'], k_aug);
end

conv = false;
inner = p;
for outiter = 1:outer

    u1 = U(:, k+1);
    v1 = V(:, k+1);
    Mu1 = MU(:, k+1);
    Nv1 = NV(:, k+1);

    if AisFunc
        Mu = A(v1, 1);
        Nv = A(u1, 2);
    else
        Mu = A * v1;
        Nv = A' * u1;
    end

    dk1 = 1;
    dk2 = -1;
    deltak1 = 0;
    gxk1 = u1;
    gyk1 = zeros(size(v1));
    gxk2 = zeros(size(u1));
    gyk2 = v1;

    for j = 1:k
        u = U(:, j);
        v = V(:, j);

        gamma2 = T(j, k+1);
        beta2 = T(k+1, j);

        alpha = T(j, j);

        d1 = 1;
        delta1 = alpha' / d1;
        d2 = -1 - abs(delta1)^2 * d1;

        eta2 = gamma2 / d1;
        sigma2 = beta2 / d2;
        lam2 = -(d1 * delta1' * eta2) / d2;

        dk1 = dk1 - abs(sigma2)^2 * d2;
        deltak1 = deltak1 - d2 * lam2 * sigma2';
        dk2 = dk2 - abs(eta2)^2 * d1 - abs(lam2)^2 * d2;

        gxk1 = gxk1 + (sigma2 * delta1)' * u;
        gyk1 = gyk1 - sigma2' * v;

        gxk2 = gxk2 - (eta2' - (lam2 * delta1)') * u;
        gyk2 = gyk2 - lam2' * v;
    end

    Mu = Mu - MU(:, 1:k) * T(1:k, k+1);
    Nv = Nv - NV(:, 1:k) * T(k+1, 1:k)';

    alpha = u1' * Mu;
    Mu = Mu - alpha * Mu1;
    Nv = Nv - alpha' * Nv1;

    deltak1 = (alpha' + deltak1) / dk1;
    dk2 = dk2 - abs(deltak1)^2 * dk1;

    pik1 = beta1 / dk1;
    pik2 = (gamma1 - deltak1 * beta1) / dk2;

    gxk2 = gxk2 - deltak1' * gxk1;
    gyk2 = gyk2 - deltak1' * gyk1;

    x = x + pik1 * gxk1 + pik2 * gxk2;
    y = y + pik1 * gyk1 + pik2 * gyk2;

    if MisFunc; u = M(Mu, 2); else; u = M \ Mu; end
    if NisFunc; v = N(Nv, 2); else; v = N \ Nv; end

    beta1 = sqrt(u' * Mu);
    gamma1 = sqrt(v' * Nv);

    u1 = u / beta1;
    v1 = v / gamma1;

    Mu0 = Mu1;
    Nv0 = Nv1;
    Mu1 = Mu / beta1;
    Nv1 = Nv / gamma1;

    T(k+1, k+1) = alpha;

    gx_1 = gxk1;
    gx0 = gxk2;
    gy_1 = gyk1;
    gy0 = gyk2;

    d_1 = dk1;
    d0 = dk2;
    delta0 = deltak1;
    pi_1 = pik1;
    pi0 = pik2;

    xi1 = pi_1 - delta0' * pi0;
    xi2 = pi0;

    j = k + 2;
    while j <= inner

        if exres
            nr = exact_norm(A, M, N, b, c, x, y);
        else
            nr = hypot(beta1 * xi2, gamma1 * xi1);
        end

        idx_res = idx_res + 1;
        resvec(idx_res) = nr;

        if nr <= tol
            exitflag = 0;
            break;
        end

        if AisFunc
            Mu = A(v1, 1);
            Nv = A(u1, 2);
        else
            Mu = A * v1;
            Nv = A' * u1;
        end

        Mu = Mu - gamma1 * Mu0;
        Nv = Nv - beta1 * Nv0;

        alpha = u1' * Mu;

        Mu = Mu - alpha * Mu1;
        Nv = Nv - alpha' * Nv1;

        if conv && mreorth > 0
            Mu = Mu - MU(:, 1:mreorth) * (U(:, 1:mreorth)' * Mu);
            Nv = Nv - NV(:, 1:mreorth) * (V(:, 1:mreorth)' * Nv);
        end
        
        if ~conv
            V(:, j) = v1;
            U(:, j) = u1;
            NV(:, j) = Nv1;
            MU(:, j) = Mu1;

            T(j, j) = alpha;
            T(j, j-1) = beta1;
            T(j-1, j) = gamma1;
        end

        if ~conv && reorth_type == 2
            Mu = Mu - MU(:, 1:j) * (U(:, 1:j)' * Mu);
            Nv = Nv - NV(:, 1:j) * (V(:, 1:j)' * Nv);
        end

        if MisFunc; u = M(Mu, 2); else; u = M \ Mu; end
        if NisFunc; v = N(Nv, 2); else; v = N \ Nv; end

        beta2 = sqrt(u' * Mu);
        gamma2 = sqrt(v' * Nv);

        sig1 = beta1 / d0;
        eta1 = gamma1 / d_1;
        lam1 = -(eta1 * delta0' * d_1) / d0;
        d1 = 1 - abs(sig1)^2 * d0;
        delta1 = (alpha' - lam1 * sig1' * d0) / d1;
        d2 = -1 - abs(eta1)^2 * d_1 - abs(lam1)^2 * d0 - abs(delta1)^2 * d1;

        pi1 = -(sig1 * pi0 * d0) / d1;
        pi2 = -(delta1 * pi1 * d1 + lam1 * pi0 * d0 + eta1 * pi_1 * d_1) / d2;

        gx1 = u1 - sig1' * gx0;
        gx2 = -delta1' * gx1 - lam1' * gx0 - eta1' * gx_1;

        gy1 = -sig1' * gy0;
        gy2 = v1 - delta1' * gy1 - lam1' * gy0 - eta1' * gy_1;

        x = x + pi1 * gx1 + pi2 * gx2;
        y = y + pi1 * gy1 + pi2 * gy2;

        xi1 = pi1 - delta1' * pi2;
        xi2 = pi2;

        % update for next iteration
        u1 = u / beta2;
        v1 = v / gamma2;
        Mu0 = Mu1;
        Nv0 = Nv1;
        Mu1 = Mu / beta2;
        Nv1 = Nv / gamma2;

        gamma1 = gamma2;
        beta1 = beta2;

        d_1 = d1;
        d0 = d2;
        pi_1 = pi1;
        pi0 = pi2;
        delta0 = delta1;

        gx_1 = gx1;
        gx0 = gx2;
        gy_1 = gy1;
        gy0 = gy2;

        j = j + 1;
    end

    if exitflag == 0 || exitflag == 2
        break;
    end

    k = k_aug;

    [U_T, S_T, V_T] = svd(T);
    Smax = S_T(1, 1);

    I = [1:k1, p-k2+1:p];
    U_T = U_T(:, I);
    V_T = V_T(:, I);
    S_T = S_T(I, I);

    T(1:k, 1:k) = S_T(1:k, 1:k);
    U(:, 1:k+1) = [U * U_T, u1];
    V(:, 1:k+1) = [V * V_T, v1];
    MU(:, 1:k+1) = [MU * U_T, Mu1];
    NV(:, 1:k+1) = [NV * V_T, Nv1];

    res1 = beta1 * V_T(p, :);
    res2 = gamma1 * U_T(p, :);

    T(1:k, k+1) = res2';
    T(k+1, 1:k) = res1;

    % [conv, absres, relres] = checkConv(abs(res1), abs(res2), k_aug, stol, Smax);
    I = (abs(res1) <= stol);
    len_conv = length( find(abs(res2(I)) <= stol) );

    if len_conv == k
        if show_info
            absres = max(abs(res1(k)), abs(res2(k)));
            fprintf(['singular triplets converge at %dth cycle, ' ...
                'absolute residual = %e, and relative residual = %e.\n'], ...
                outiter, absres, absres/Smax);
        end

        conv = true;
        exitflag = 2;
        inner = maxit;
    end

    beta1 = -beta1 * xi2;
    gamma1 = -gamma1 * xi1;
end

if show_info
    if exitflag == 1
        fprintf(['Reach maximum cycle: %d. No covergence to required ' ...
            'tolerance %e. Absolute error is %.4e.\n'], outiter, tol, resvec(idx_res));
    end

    if exitflag == 2
        warning(['The target elliptic singular triplets converges, ' ...
            'reaching the maximum iterations at the non-restarting stage. ' ...
            'No covergence to required tolerance %e. Absolute error is %.4e.\n'], ...
            tol, resvec(idx_res));
    end
end

resvec = resvec(1:idx_res);

[varargout{1:5}] = deal(U(:, 1:k+1), T(1:k+1, 1:k+1), V(:, 1:k+1), MU(:, 1:k+1), NV(:, 1:k+1));

end