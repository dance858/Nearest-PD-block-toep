function [lambda, lb_opt, U1, U2, num__fstchol_called] = initial_guess(C, m, p, lambda0)
% IN:
%      C: matrix of size m x (p+1)m
%      lambda0: initial guess for a value lambda such that Toep(C^T) + lambda*I
%               is positive definite.
% OUT:
%      lambda: value such that that Toep(C^T) + lambda*I is positive
%              definite.
%      lb_opt: lower bound on the optimal value of lambda.
%      L1 and L2 are upper triangular matrices satisfying that
%      Toep(C^T) + lambda*I = L1*L1' and (Toep(C^T) + lambda*I)^{-1} = L2*L2'.

% Attempt to do a Cholesky factorization for the initial guess. 
pd_flag = true;
C(:, 1:m) = C(:, 1:m) + lambda0*eye(m);
try
    [U1,~,U2] = fstchol(C');
catch
   pd_flag = false; 
end
C(:, 1:m) = C(:, 1:m) - lambda0*eye(m);

% If the Cholesky factorization succeeds, then we can initialize at 
% lambda0. If the Cholesky factorization fails, then we initialize
% using Gershgorins disc thm. 
if pd_flag
    lb_opt = -min(diag(C(:, 1:m)));
    lambda = lambda0;
    num__fstchol_called = 1;
else
    lb_opt = lambda0;
    % Gershgorin
    TC = block_toep(mat2cell(C', m*ones(1, p+1), m)');
    col_sums = sum(abs(TC), 1);
    diagonal_elements = diag(TC)+ abs(diag(TC));
    lambda = -min(diagonal_elements' - col_sums);
    C(1:m, 1:m) = C(1:m, 1:m) + lambda*eye(m);
    num__fstchol_called = 2;
    try
        [U1,~,U2] = fstchol(C');
    catch
        error("Cholesky in initial_guess failed. This should never happen. ") 
    end
end

 if lb_opt > lambda
      error("Inconsistent initial guess. This should never happen. "); 
 end
end