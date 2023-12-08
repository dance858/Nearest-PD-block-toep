clear; clc;
addpath('../algorithms/')
addpath('../algorithms/entropic_prox_operator')

%% Settings for experiment
m = 5;
all_p = [20, 40, 60, 80, 100];
num_runs = 5;

%% Settings for IGA
max_iter = 1000; tol_subp = 1e-8; L0 = 5;
beta = 2; gamma = 1.05; c = 0.001;
verbose = true;

for p = all_p
    for run = 2:num_runs
        filename = "data/m=" + num2str(m) + "/BT_m=" + num2str(m) + "_p=" + ...
                    num2str(p) + "_run=" + num2str(run) + ".mat";
        load(filename, "BT");

        % Extract first block column as row
        A = zeros(m, m*(p+1));
        A(:, 1:m) = BT(1:m, 1:m);
        for k = 0:p
            A(:, k*m + 1:(k+1)*m) = BT(k*m + 1:(k+1)*m, 1:m);
        end
        
        % Solve with IGA
        [Z_IGA, f_conj_val_Z, stats_proj, all_L, objs, opt_lambdas, time,  ...
         UB_primal, LB_primal, time_MIS_proj, time_iter, total_num_of_chol_each_iter] = ...
        solve_dual_IGA(A, max_iter, tol_subp, L0, beta, gamma, c, verbose);
    
        fprintf("Number of iterations: %i \n", length(objs));
        fprintf("Average number of Newton steps: %f \n", mean(stats_proj(1, :)));
        fprintf("Max number of Cholesky per prox: %f: \n", max(stats_proj(2, :)));
        fprintf("Max number of Cholesky per iter: %f: \n \n  ", max(total_num_of_chol_each_iter));
    
        filename = "data_IGA/m=" + num2str(m) + "/BT_m=" + num2str(m) + "_p=" + num2str(p) + "_run=" + ...
                    num2str(run) + "_c=" + num2str(c) + "_L=" + num2str(L0) + ".mat";
        save(filename);
        
    end
end


