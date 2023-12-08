function [Z, f_conj_val_Z, stats_proj, all_L, objs, opt_lambdas, time,  ...
         UB_primal, LB_primal, time_MIS_proj, time_iter, total_num_of_chol_each_iter] = ...
            solve_dual_IGA(A, max_iter, tol_subp, L0, beta, gamma, c, verbose)
% 
% IN:
%     A: The first block column of the block Toeplitz matrix. 
%        Size m x m(p+1).

m = size(A, 1); p = size(A, 2)/m - 1;
Z = zeros(m, (p+1)*m);
Z(:, 1:m) = eye(m);
Z_val_phi = 0 + c*m^2;
Z_grad_phi = zeros(m, m*(p+1));                 
Z_grad_phi(:, 1:m) = -eye(m) + 2*c*m*eye(m); 

[Z, f_conj_val_Z, stats_proj, all_L, objs, opt_lambdas, time, UB_primal, LB_primal, ...
    time_MIS_proj, time_iter, total_num_of_chol_each_iter] = ...
   IGA(Z, Z_grad_phi, Z_val_phi, A, max_iter, tol_subp, L0, beta, gamma, c, verbose);
     
end