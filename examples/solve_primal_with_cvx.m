function [TX, cvx_optval, solve_time] = solve_primal_with_cvx(A_CVX, m, p, solver)
% Computes TX - the positive semidefinite block Toeplitz matrix
%               that is nearest to T(A), which is the block Toeplitz matrix
%               whose first column is A_CVX.
% IN:
%        A_CVX: Represents the first block column of T(A). 
%               Size m(p+1) x m. 
% 

block_toep_A = block_toep(mat2cell(A_CVX, m*ones(1, p+1), m)');

if solver == "scs"
   cvx_solver scs
else
   cvx_solver SDPT3  
end
cvx_begin sdp
    variable X_CVX(m*(p+1), m)      % first block column
    variable T(m*(p+1), m*(p+1)) symmetric; 
    obj = square_pos(norm(T - block_toep_A, 'fro'));
    minimize(obj);
    subject to
       T >= 0;
       for row = 0:p
           for col = row:p
               T(row*m + 1:(row+1)*m, col*m+1:(col+1)*m) == ...
                   X_CVX((col - row)*m + 1:(col - row + 1)*m, :)';
           end
       end
start = tic;
cvx_end
solve_time = toc(start);

TX = block_toep(mat2cell(X_CVX, m*ones(1, p+1), m)');
end