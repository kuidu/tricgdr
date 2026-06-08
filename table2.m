clc; close all; clearvars;

mats = {'gupta3', 'g7jac060sc', 'rajat27', 'TSOPF_RS_b300_c2'};
% ks = [10, 15, 20, 25];
ks = [10, 20, 40, 40];
ps = [80, 80, 100, 120];

sz = [length(mats), length(ks)+2];
var_types = cell(1, sz(2));
[var_types{:}] = deal('double');
var_names = ["(p="+string(ps)+",k="+string(ks)+")", "TriCG", "MINRES"];
sum_data = table('Size', sz, 'VariableTypes', var_types, ...
                 'VariableNames', var_names, 'RowNames', mats');

tol = 1e-8;
stol = 1e-10;
maxcycle = 10;
maxit = 80000;

% run 'it' times, and take the average time
it = 10;
%%
% the average runtime of TriCG
for i = 1:length(mats)
    mat_name = mats{i};
    load(['mats', filesep, mat_name, '.mat']);
    A = Problem.A;

    [m, n] = size(A);

    b = ones(m ,1);
    c = ones(n, 1);
    b = b / norm(b);
    c = c / norm(c);

    pause(300)

    tricg_time = tic;
    for j = 1:it
        [~, ~, flag, ~] = tricg(A, b, c, [], [], tol, maxit);
        if flag ~= 0
            fprintf('tricg fails to converge for the matrix: %s\n', mat_name);
        end
    end
    t1 = toc(tricg_time);
    t1 = t1 / it;
    sum_data{i, end-1} = round(t1, 2);

end
disp(rows2vars(sum_data(:, end-1))); %[output:2b19ef5e]
%%
% the average time of MINRES
for i = 1:length(mats)
    mat_name = mats{i};
    load(['mats', filesep, mat_name, '.mat']);
    A = Problem.A;

    [m, n] = size(A);

    b = ones(m ,1);
    c = ones(n, 1);
    b = b / norm(b);
    c = c / norm(c);

    K = [speye(m), A; A', -speye(n)];
    f = [b; c];

    pause(300)

    minres_time = tic;
    for j = 1:it
        [~, ~, resvec3] = minres(K, f, [], tol, maxit);
        if flag ~= 0
            fprintf('minres fails to converge for the matrix: %s\n', mat_name);
        end
    end
    t2 = toc(minres_time);
    t2 = t2 / it;
    sum_data{i, end} = round(t2, 2);
end
disp(rows2vars(sum_data(:, end))); %[output:2e52c5ef]
%%
% the average runtime of TriCG-DR

for i = 1:length(mats)
    mat_name = mats{i};
    load(['mats', filesep, mat_name, '.mat']);
    A = Problem.A;

    [m, n] = size(A);

    b = ones(m ,1);
    c = ones(n, 1);
    b = b / norm(b);
    c = c / norm(c);

    K = [speye(m), A; A', -speye(n)];
    f = [b; c];

    for j = 1:length(ks)
        k = ks(j);
        p = ps(j);
        
        pause(300)

        tricgdr_time = tic;
        for jj = 1:it
            [~, ~, flag, ~] = tricgdr(A, b, c, [], [], k, p, tol, maxcycle, ...
                maxit=maxit, stol=stol);
            if flag ~= 0
                fprintf('tricg-dr fails to converge for the matrix: %s, k: %d\n', mat_name, k);
            end
        end
        t3 = toc(tricgdr_time);
        t3 = t3 / it;
        sum_data{i, j} = round(t3, 2);
    end

end
disp(sum_data(1:end, 1:end-2)); %[output:5bc9e2d5]
%%
writetable(sum_data, 'table2.csv', Delimiter=' ', WriteRowNames=true);
disp(sum_data); %[output:5edeca42]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:2b19ef5e]
%   data: {"dataType":"text","outputData":{"text":"    <strong>OriginalVariableNames<\/strong>    <strong>gupta3<\/strong>    <strong>g7jac060sc<\/strong>    <strong>rajat27<\/strong>    <strong>TSOPF_RS_b300_c2<\/strong>\n    <strong>_____________________<\/strong>    <strong>______<\/strong>    <strong>__________<\/strong>    <strong>_______<\/strong>    <strong>________________<\/strong>\n\n          {'TriCG'}          17.88       11.21        18.53          23.09      \n\n","truncated":false}}
%---
%[output:2e52c5ef]
%   data: {"dataType":"text","outputData":{"text":"    <strong>OriginalVariableNames<\/strong>    <strong>gupta3<\/strong>    <strong>g7jac060sc<\/strong>    <strong>rajat27<\/strong>    <strong>TSOPF_RS_b300_c2<\/strong>\n    <strong>_____________________<\/strong>    <strong>______<\/strong>    <strong>__________<\/strong>    <strong>_______<\/strong>    <strong>________________<\/strong>\n\n         {'MINRES'}          31.86       10.05         10            30.19      \n\n","truncated":false}}
%---
%[output:5bc9e2d5]
%   data: {"dataType":"text","outputData":{"text":"                        <strong>(p=80,k=10)<\/strong>    <strong>(p=80,k=20)<\/strong>    <strong>(p=100,k=40)<\/strong>    <strong>(p=120,k=40)<\/strong>\n                        <strong>___________<\/strong>    <strong>___________<\/strong>    <strong>____________<\/strong>    <strong>____________<\/strong>\n\n    <strong>gupta3          <\/strong>       13.84          11.51           9.54            9.48    \n    <strong>g7jac060sc      <\/strong>        8.42           8.98          11.84           11.98    \n    <strong>rajat27         <\/strong>        9.75           8.43           5.88            5.71    \n    <strong>TSOPF_RS_b300_c2<\/strong>       18.24          17.02          12.39           15.72    \n\n","truncated":false}}
%---
%[output:5edeca42]
%   data: {"dataType":"text","outputData":{"text":"                        <strong>(p=80,k=10)<\/strong>    <strong>(p=80,k=20)<\/strong>    <strong>(p=100,k=40)<\/strong>    <strong>(p=120,k=40)<\/strong>    <strong>TriCG<\/strong>    <strong>MINRES<\/strong>\n                        <strong>___________<\/strong>    <strong>___________<\/strong>    <strong>____________<\/strong>    <strong>____________<\/strong>    <strong>_____<\/strong>    <strong>______<\/strong>\n\n    <strong>gupta3          <\/strong>       13.84          11.51           9.54            9.48        17.88    31.86 \n    <strong>g7jac060sc      <\/strong>        8.42           8.98          11.84           11.98        11.21    10.05 \n    <strong>rajat27         <\/strong>        9.75           8.43           5.88            5.71        18.53       10 \n    <strong>TSOPF_RS_b300_c2<\/strong>       18.24          17.02          12.39           15.72        23.09    30.19 \n\n","truncated":false}}
%---
