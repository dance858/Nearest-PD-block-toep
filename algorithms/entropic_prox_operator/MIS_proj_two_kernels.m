function [sol, grad_phi, phi_1, phi_2, stats, opt_lambda] =  ...
    MIS_proj_two_kernels(C, c, tol, lambda0, kernel)
% Let phi_1(X) = -1/(2*pi)*int_0^{2\pi} log det F_X(e^{j w}) dw
% and phi_2(X) = c*Tr(X0^2).
%
% If kernel == 1, this function solves the problem
% minimize      <C, Y> + phi_1(Y)
% subject to    Tr(Y0) = c.
% If kernel == 2, this function solves the problem
% minimize <C, Y> + phi_1(Y) + c*t^2
% subject to Tr(Y0) = t
%  Note: * Y = [Y0, Y1, ..., Yp], C = [C0, C1, ..., Cp]
%        * <C, Y> = Tr(C0*Y0) + 2*[Tr(C1'*Y1) + ... + Tr(Cp'*Yp)] etc.
%        * phi_1(X) = -1/(2*pi)*int_0^{2\pi} log det F_X(e^{j w}) dw
%        * phi_2(X) = c*Tr(X0^2)
%  INPUT:
%        C: matrix of size m x m(p+1)
%        c: used in the definition of phi_2
%      tol: tolerance for Newton's method 
%  lambda0: initial guess for dual problem.

maxiter = 30; m = size(C, 1); p = size(C, 2)/m - 1;
E = zeros(m, (p+1)*m); E(:, 1:m) = eye(m);
stats = zeros(2,1); %[num of iter, num of times fstchol is called]

% Compute initial guess for dual variable lambda, and a lower bound on 
% the optimal value of lambda. Also compute upper triangular U such that
% (Toep(C^T) + lambda*I)^{-1} = U*U'.
[lambda, lb_opt, ~, U, num__fstchol_called] = initial_guess(C, m, p, lambda0);
stats(2) = stats(2) + num__fstchol_called;

% Evaluate dual objective function. To evaluate the log-determinant in 
% a numerically safe way we use Cholesky. 
s_val = 2*sum(log(diag(chol(U(1:m, :)*U(1:m, :)'))));
if kernel == 1
    s_val = s_val + c*lambda;
else  
    s_val = s_val + lambda^2/(4*c);  
end

for ii = 1:maxiter
    if ii == maxiter - 1
       disp('Max iterations in MISproj!') 
    end
    
    if ii >= 20
       w = 1; 
    end
    
    % Compute first and second derivative of s in the current iterate.
    [s_prime, s_prime_prime] = s_diff(U, m, c, lambda, kernel);
   
    % Check termination criteria (Newton decrement)
    if abs(s_prime^2/(2*s_prime_prime)) < tol 
        break;
    end                           
     
    % Take a Newton step. 
    dir = -s_prime/s_prime_prime;
    [step_size, U, s_val, lb_opt, num__fstchol_called] = ...
        compute_step_size(lambda, dir, lb_opt, s_prime, s_val, C, E, c, kernel);
    lambda = lambda + step_size*dir;
    stats(2) = stats(2) + num__fstchol_called;
end

% First recover the spectral factor B of the optimal Y, and then Y itself. 
% Note that R has in fact already been computed inside 'compute_s_derivatives',
% so possible to refactor. 
R = U*U(1:m, :)'; 
B = zeros(m*(p+1), m);
B(1:m, :) = chol(R(1:m, :))';
B(m:end, :) = R(m:end, :)/B(1:m, :)';
sol = block_diag_sum_two_transpose(B*B', m, p);

% Recover the value and gradient of phi at the optimal Y.
% When evaluating phi we take into account that B0 is upper triangular.
opt_lambda = lambda;
grad_phi = -C;
phi_1 = -2*sum(log(diag(B(1:m, :))));
phi_2 = c*trace(sol(:, 1:m))^2;

%if kernel == 1
%   grad_phi(:, 1:m) = grad_phi(:, 1:m) - opt_lambda*eye(m); 
%   phi_1 = -2*sum(log(diag(B(1:m, :))));
%   phi_2 = -1337;
%else
%   grad_phi(:, 1:m) = grad_phi(:, 1:m) + (2*c*trace(sol(:, 1:m)) - opt_lambda)*eye(m);
%   phi_1 = -2*sum(log(diag(B(1:m, :))));
%   phi_2 = opt_lambda^2/(4*c);  
   %phi_2 = c*trace(sol(:, 1:m))^2;
%end

% Update statistics for solver.
stats(1) = ii;

end