function [UB, feasible_X] = compute_UB(A, Z, m, p)

X = A;
for k = 0:p 
  X(:, k*m + 1:(k+1)*m) = X(:, k*m + 1:(k+1)*m) + 0.5/(p+1-k)*Z(:, k*m+1:(k+1)*m);
end

% Check if primal iterate is feasible
try
    [chol_factor] = fstchol(X');
catch
   TX = block_toep(mat2cell(X, m, m*ones(1, p+1))');
   X(:, 1:m) = X(:, 1:m) - min(eig(TX))*eye(m);
end

feasible_X = X;
UB = norm(block_toep(mat2cell(X-A, m, m*ones(1, p+1))'), 'fro')^2;
end