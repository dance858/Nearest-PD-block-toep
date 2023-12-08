function [feasible, UB, X] = primal_feas(A, Z, m, p)

X = A;
for k = 0:p 
  X(:, k*m + 1:(k+1)*m) = X(:, k*m + 1:(k+1)*m) + 0.5/(p+1-k)*Z(:, k*m+1:(k+1)*m);
end

% X_cell = {};
% for k = 0:p
% X_cell{k+1} =  X(:, k*m + 1:(k+1)*m);
% end
% TX = block_toep(X_cell);

% Check if primal iterate is feasible
try
    [chol_factor] = fstchol(X');
    feasible = true;
    UB = 1/(p+1)*norm(Z(:, 1:m), 'fro')^2;
    for k = 1:p
       UB = UB + 2/(p + 1 - k)*norm(Z(:, m*k+1:(m+1)*k), 'fro')^2; 
    end
    UB = UB/4;
catch
   feasible = false;
   UB = 1e-5;
   X_cell = {};
   for k = 0:p
      X_cell{k+1} =  X(:, k*m + 1:(k+1)*m);
   end
   TX = block_toep(X_cell);
   fprintf("Smallest eigenvalue TX: %f \n", min(eig(TX)));
end


end