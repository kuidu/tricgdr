clc; close all; clearvars;

sz = [1, 4];
var_types = cell(1, sz(2));
[var_types{:}] = deal('double');
var_names = ["TriCG", "MINRES", "TriCG-DR+D-TriCG", "TriCG-DR"];
sum_data = table('Size', sz, 'VariableTypes', var_types, ...
    'VariableNames', var_names, 'RowNames', "CPU time");
%%
load('mats/ifiss.mat');

tau = 1e-1;
M = Mh - tau * Ah;
C = C - 1e-10*speye(size(C));
N = -tau * C;
A = -tau * B;

b0 = tau * f;
c0 = -tau * g;

dM = decomposition(M, 'chol');
dN = decomposition(N, 'chol');
Mfunc = @(x, t) Mopertor(M, dM, x, t);
Nfunc = @(x, t) Mopertor(N, dN, x, t);

tol = 1e-10;
maxit = 4000;

stol = 1e-10;
k = 40;
p = 120;
maxcycle = 10;
maxaug = 4000;

num_rhs = 10;
it = 10;
%%
% TriCG for all the systems
pause(300)
tricg_time = tic;
for l = 1:it
    [x, y] = tricg(A, b0, c0, Mfunc, Nfunc, tol, maxit);
    for i = 2:num_rhs
        b = b0 + Mh * x;
        [x, y, flag] = tricg(A, b, c0, Mfunc, Nfunc, tol, maxit);
        
        if flag ~= 0
            fprintf('TriCG fails to converge for the %dth system. \n', i)
        end
    end
end
t1 = toc(tricg_time);
t1 = t1 / it;
sum_data{1, 1} = round(t1, 2);
disp(sum_data); %[output:0b87ee6f]
%%
% MINRES for all the systems
K = [M, A; A', -N];
m = size(A, 1);
pre = @(x, t) [Mfunc(x(1:m), t); Nfunc(x(m+1:end), t)];

pause(300)
minres_time = tic;
for l = 1:it
    x = minres(K, [b0; c0], pre, tol, maxit);
    for i = 2:num_rhs
        b = b0 + Mh * x(1:m);
        [x, flag] = minres(K, [b; c0], pre, tol, maxit);

        if flag ~= 0
            fprintf('MINRES fails to converge for the %dth system. \n', i)
        end
    end
end
t2 = toc(minres_time);
t2 = t2 / it;
sum_data{1, 2} = round(t2, 2);
disp(sum_data); %[output:20f87ee7]
%%
pause(300)
dtricg_time = tic;
for l = 1:it
    % TriCG-DR for the 1st system
    [x, y, ~, ~, U, T, V, MU, NV] = tricgdr(A, b0, c0, Mfunc, Nfunc, k, p, tol, ...
        maxcycle, maxit=maxit, stol=stol);

    % D-TriCG for the 2nd to the last systems
    for i = 2:num_rhs
        b = b0 + Mh * x;
        d1 = [eye(k), T(1:k, 1:k); T(1:k, 1:k)', -eye(k)] \ [U(:,1:k)'*b; V(:,1:k)'*c0];
        d2 = [eye(k+1, k), T(1:k+1, 1:k); T(1:k, 1:k+1)', -eye(k+1, k)] * d1;
        rx = b - MU(:, 1:k+1) * d2(1:k+1);
        ry = c0 - NV(:, 1:k+1) * d2(k+2:end);

        [dx, dy, flag] = tricg(A, rx, ry, Mfunc, Nfunc, tol, maxit, augU=U(:,1:k), ...
            augMU=MU(:,1:k), augV=V(:,1:k), augNV=NV(:,1:k));
        x = U(:, 1:k) * d1(1:k) + dx;
        y = V(:, 1:k) * d1(k+1:end) + dy;

        if flag ~= 0
            fprintf('D-TriCG fails to converge for the %dth system. \n', i)
        end
    end
end
t3 = toc(dtricg_time);
t3 = t3 / it;
sum_data{1, 3} = round(t3, 2);
disp(sum_data); %[output:3cf04b51]
%%
% TriCG-DR for all the systems
pause(300)
tricgdr_time = tic;
for l = 1:it
    [x, y] = tricgdr(A, b0, c0, Mfunc, Nfunc, k, p, tol, ...
        maxcycle, maxit=maxit, stol=stol);
    for i = 2:num_rhs
        b = b0 + Mh * x;
        [x, y, flag] = tricgdr(A, b, c0, Mfunc, Nfunc, k, p, tol, ...
            maxcycle, maxit=maxit, stol=stol);

        if flag ~= 0
            fprintf('TriCG-DR fails to converge for the %dth system. \n', i)
        end
    end
end
t4 = toc(tricgdr_time);
t4 = t4 / it;
sum_data{1, 4} = round(t4, 2);
disp(sum_data); %[output:608f79b9]
%%
writetable(sum_data, 'table3.csv', Delimiter=' ', WriteRowNames=true);
disp(sum_data); %[output:15c2e8de]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:0b87ee6f]
%   data: {"dataType":"text","outputData":{"text":"                <strong>TriCG<\/strong>     <strong>MINRES<\/strong>    <strong>TriCG-DR+D-TriCG<\/strong>    <strong>TriCG-DR<\/strong>\n                <strong>______<\/strong>    <strong>______<\/strong>    <strong>________________<\/strong>    <strong>________<\/strong>\n\n    <strong>CPU time<\/strong>    144.43      0              0               0    \n\n","truncated":false}}
%---
%[output:20f87ee7]
%   data: {"dataType":"text","outputData":{"text":"                <strong>TriCG<\/strong>     <strong>MINRES<\/strong>    <strong>TriCG-DR+D-TriCG<\/strong>    <strong>TriCG-DR<\/strong>\n                <strong>______<\/strong>    <strong>______<\/strong>    <strong>________________<\/strong>    <strong>________<\/strong>\n\n    <strong>CPU time<\/strong>    144.43    217.25           0               0    \n\n","truncated":false}}
%---
%[output:3cf04b51]
%   data: {"dataType":"text","outputData":{"text":"                <strong>TriCG<\/strong>     <strong>MINRES<\/strong>    <strong>TriCG-DR+D-TriCG<\/strong>    <strong>TriCG-DR<\/strong>\n                <strong>______<\/strong>    <strong>______<\/strong>    <strong>________________<\/strong>    <strong>________<\/strong>\n\n    <strong>CPU time<\/strong>    144.43    217.25         78.85             0    \n\n","truncated":false}}
%---
%[output:608f79b9]
%   data: {"dataType":"text","outputData":{"text":"                <strong>TriCG<\/strong>     <strong>MINRES<\/strong>    <strong>TriCG-DR+D-TriCG<\/strong>    <strong>TriCG-DR<\/strong>\n                <strong>______<\/strong>    <strong>______<\/strong>    <strong>________________<\/strong>    <strong>________<\/strong>\n\n    <strong>CPU time<\/strong>    144.43    217.25         78.85           121.41 \n\n","truncated":false}}
%---
%[output:15c2e8de]
%   data: {"dataType":"text","outputData":{"text":"                <strong>TriCG<\/strong>     <strong>MINRES<\/strong>    <strong>TriCG-DR+D-TriCG<\/strong>    <strong>TriCG-DR<\/strong>\n                <strong>______<\/strong>    <strong>______<\/strong>    <strong>________________<\/strong>    <strong>________<\/strong>\n\n    <strong>CPU time<\/strong>    144.43    217.25         78.85           121.41 \n\n","truncated":false}}
%---
