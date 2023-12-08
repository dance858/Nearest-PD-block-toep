function [Z, f_conj_val_Z, stats_proj, all_L, objs, opt_lambdas, time, UB_primal, ...
         LB_primal, time_MIS_proj, time_iter, total_num_of_chol_each_iter] = ...
         IGA(Z, Z_grad_phi, Z_val_phi, A, max_iter, tol_subp, L0, beta, gamma, c, verbose)
     
% An improved interior gradient method for solving
%      minimize     f*(Z)
%      subject to   Z \in K*
% where 
%     f*(Z) = 1/4*||C^{1/2} Z||^2 + <A, Z> 
%     K* = cone of nonnegative trigonometric matrix polynomials
%
% INPUT:
%               Z:  Starting point, size m x m(p+1)
%      Z_grad_phi:  Gradient of kernel. Note that it should be full kernel 
%                   (so phi_1 + phi_2).
%      Z_val_phi:   Value of (full) kernel.
%              A:   Problem data.
%       max_iter: the algorithm terminates if max_iter is reached      
%       tol_subp: the tolerance used for evaluating the entropic
%                 prox operator. Typical value is 1e-7
%        max_iter: the algorithm terminates if max_iter is reached 
%              L0: inital estimate of smoothness. Should be large?
%     beta, gamma: line search parameters > 1.

tic;
% Initialization
m = size(Z, 1); p = size(Z, 2)/m - 1; kernel = 2;
V = Z; V_grad_phi = Z_grad_phi; V_phi = Z_val_phi;
f_conj_val_Z = val_f_conj(Z, A, m, p);

% beta = 1.5; gamma = 1.1;
opt_lambda = 0;   % Initial guess for projection problem.
L = L0*gamma;     % L corresponds to L_k.
L_old = -1;       % L_old corresponds to L_{k-1}      
theta_old = -1;   % theta_old corresponds to theta_{k-1}

% Quantities for tracking the progress.
stats_proj = zeros(2, 0); all_L = [];
objs = [f_conj_val_Z]; opt_lambdas = []; 
LB_primal = - f_conj_val_Z;
[UB_primal, ~] = compute_UB(A, Z, m, p);
time_MIS_proj = [];
time_iter = [];
total_num_of_chol_each_iter = [];

for k = 1:max_iter
    start = tic;
    total_num_of_chol_each_iter(end+1) = 0;
    if verbose && ~rem(k, 10)
        fprintf('iter/L/obj: %i / %f / %f \n ', k, L, objs(end))
    end
    %fprintf('iter/L/obj: %i / %f / %f \n ', k, L, objs(end))
    
    L = L/gamma;              
    terminate_linesearch = false;
    
    while ~terminate_linesearch
       if k == 1
          theta = 1; 
       else
          theta = -theta_old^2*L_old/(2*L) +...
                  L_old/2*sqrt(theta_old^4/L^2 + 4*theta_old^2/(L*L_old));
       end
       tau = 1/(L*theta);
       Y = (1-theta)*Z + theta*V;   
       f_conj_grad_Y = grad_f_conj(Y, A, m, p);   % gradient in Y
       
       MIS_proj_start = tic;
       [V_cand, V_grad_phi_cand, phi_1_V_cand, phi_2_V_cand, stats, opt_lambda] = ...
           MIS_proj_two_kernels(tau*f_conj_grad_Y - V_grad_phi, c, tol_subp, ...
                                opt_lambda, kernel); 
       time_MIS_proj(end+1) = toc(MIS_proj_start);
       total_num_of_chol_each_iter(k) = total_num_of_chol_each_iter(k) + stats(2);
                            
       Z_cand = (1-theta)*Z + theta*V_cand;
       
       % Store some statistics
       stats_proj(:, end+1) = stats;
       opt_lambdas(end+1) = opt_lambda;
       
       % Check descent criteria.
       bregman_dist_V_cand_and_V = (phi_1_V_cand + phi_2_V_cand) - V_phi - ...
                                   inner_prod_E(V_grad_phi, V_cand - V, m); 
       if bregman_dist_V_cand_and_V < 0
          error("Negative Bregman divergence") 
       end
       rhs_descent_test = (1-theta)*f_conj_val_Z ...
            + theta*(val_f_conj(Y, A, m, p) + inner_prod_E(f_conj_grad_Y, V_cand - Y, m)...
                     + 1/tau*bregman_dist_V_cand_and_V);
       
       f_conj_val_Z_cand = val_f_conj(Z_cand, A, m, p);           
       
       %if true
       if f_conj_val_Z_cand >= rhs_descent_test
          q = 1; 
       end
       
       if f_conj_val_Z_cand <= rhs_descent_test 
            terminate_linesearch = true;
            L = L/beta;     % To compensate for later
       end
       L = L*beta;
    end
    Z = Z_cand;
    V = V_cand;
    V_phi = phi_1_V_cand + phi_2_V_cand;
    V_grad_phi = V_grad_phi_cand;
    f_conj_val_Z = f_conj_val_Z_cand;
    L_old = L;    
    theta_old = theta;
    objs(end+1) = f_conj_val_Z;
    LB_primal = max(LB_primal, -f_conj_val_Z);
    
    all_L(end+1) = L;
    time_iter(end+1) = toc(start);
    % TODO: Check termination criteria
    if ~rem(k, 20)
        [UB, ~] = compute_UB(A, Z, m, p);
        UB_primal = min(UB, UB_primal);
        
        if (UB_primal - LB_primal)/UB_primal < 1e-3
            break;
        end
    end
    
end
time = toc;
end