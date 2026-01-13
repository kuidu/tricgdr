function [Uk, Sk, Vk, MUk, NVk, Tk, len_conv] = gssydr(A, b, c, M, N, k, p, tol, maxcycle, opt)
% GSSYDR will find k approximate elliptic singular triplets of A, using the
% generalized Saunders--Simon--Yip tridiagonalization process.
%
% Syntaxes
% --------------
% [U, S, V, MU, NV, T, l] = gssydr(A, b, c)
% [U, S, V, MU, NV, T, l] = gssydr(A, b, c, M, N)
% [U, S, V, MU, NV, T, l] = gssydr(A, b, c, M, N, k, p)
% [U, S, V, MU, NV, T, l] = gssydr(A, b, c, M, N, k, p, tol)
% [U, S, V, MU, NV, T, l] = gssydr(A, b, c, M, N, k, p, tol, maxcycle)
% [U, S, V, MU, NV, T, l] = gssydr(A, b, c, M, N, k, p, tol, maxcycle)
% [U, S, V, MU, NV, T, l] = gssydr(__, Name, Value)
%
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
% maxcycle      number of maximum cycles
%
%
% Name, Value   name-Value pairs determine other options. 
%                   `reorth`:   the type of reorthogonalization. `1` means
%                               only reorthogonalize one side. `2` means 
%                               reorthogonalize both two sides
%
%
% Returns
% --------------
%
% U, S, V, MU, NV, T        approximate singular triplets
%
% l                         number of converged approximate singular triplets

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
    opt.reorth      (1, 1)  double  {mustBeMember(opt.reorth, [0, 2])}      =  2
end

reorth_type = opt.reorth;

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

lb = length(b);
lc = length(c);

U = zeros(lb, p);
V = zeros(lc, p);
MU = zeros(lb, p);
NV = zeros(lc, p);
T = zeros(p, p);

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

k = 0;
k_aug = k1 + k2;

if k_aug >= p
    k_aug = min(max(p-10, 10), floor(p/2));
    warning(['dimension of augment space should be less than dimension of subspace. ' ...
        'set k = %d\n'], k_aug);
end

inner = p;
ttol = eps^(2/3); %tolerance for terminate
for outiter = 1:maxcycle

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

    Mu = Mu - MU(:, 1:k) * T(1:k, k+1);
    Nv = Nv - NV(:, 1:k) * T(k+1, 1:k)';

    alpha = u1' * Mu;
    Mu = Mu - alpha * Mu1;
    Nv = Nv - alpha' * Nv1; 

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

    j = k + 2;
    while j <= inner

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

        
        V(:, j) = v1;
        U(:, j) = u1;
        NV(:, j) = Nv1;
        MU(:, j) = Mu1;

        T(j, j) = alpha;
        T(j, j-1) = beta1;
        T(j-1, j) = gamma1;


        if reorth_type == 2
            Mu = Mu - MU(:, 1:j) * (U(:, 1:j)' * Mu);
            Nv = Nv - NV(:, 1:j) * (V(:, 1:j)' * Nv);
        end

        if MisFunc; u = M(Mu, 2); else; u = M \ Mu; end
        if NisFunc; v = N(Nv, 2); else; v = N \ Nv; end

        beta2 = sqrt(u' * Mu);
        gamma2 = sqrt(v' * Nv);

        if min(abs(beta2), abs(gamma2)) <= ttol
            p = j;
            k_aug = min(k_aug, p);
            k1 = k_aug;
            k2 = 0;
            T = T(1:j, 1:j);
            warning(['beta or gamma is approximate zero, gSSY-DR terminates, ' ...
                'and we only compute the first %d singular triplets'], k_aug);
            break;
        end

        % update for next iteration
        u1 = u / beta2;
        v1 = v / gamma2;
        Mu0 = Mu1;
        Nv0 = Nv1;
        Mu1 = Mu / beta2;
        Nv1 = Nv / gamma2;

        gamma1 = gamma2;
        beta1 = beta2;

        j = j + 1;
    end

    k = k_aug;

    [U_T, S_T, V_T] = svd(T);

    I = [1:k1, p-k2+1:p];

    U_T = U_T(:, I);
    V_T = V_T(:, I);
    S_T = S_T(I, I);

    T(1:k, 1:k) = S_T;
    U(:, 1:k+1) = [U * U_T, u1];
    V(:, 1:k+1) = [V * V_T, v1];
    MU(:, 1:k+1) = [MU * U_T, Mu1];
    NV(:, 1:k+1) = [NV * V_T, Nv1];

    res1 = beta1 * V_T(p, :);
    res2 = gamma1 * U_T(p, :);

    T(1:k, k+1) = res2';
    T(k+1, 1:k) = res1;

    % [conv, ~, ~, len_conv] = checkConv(abs(res1), abs(res2), k, tol, Smax);
    I = (abs(res1) <= tol);
    len_conv = length( find(abs(res2(I)) <= tol) );

    if len_conv == k
        % fprintf(['The %d singular triplets converge after %d cycles, ' ...
        %     'and the absoulte error is %e\n'], k, outiter, max(abs([res1, res2])) );
        break;
    end
end

if len_conv == 0
    warning(['Reaching the maximum number of cycles %d, ' ...
        'no elliptic singular triplets converge.\n'], maxcycle);
end

if len_conv > 0 && len_conv < k
    warning(['Reaching the maximum number of cycles %d, ' ...
        'only the first %d elliptic singular triplets converge, error=%e.\n'], ...
        maxcycle, len_conv, max(abs([res1, res2])) );
end

Uk = U(:, 1:k+1);
Vk = V(:, 1:k+1);
MUk = MU(:, 1:k+1);
NVk = NV(:, 1:k+1);
Sk = S_T;
Tk = T(1:k+1, 1:k+1);

end
