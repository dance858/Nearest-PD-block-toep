function [obj] = val_f_conj(Z, A, m, p)

%term2 = 2*trace(A(:, m+1:end)'*Z(:, m+1:end)) + trace(A(:, 1:m)*Z(:, 1:m));
term2 = 2*sum(A(:, m+1:end).*Z(:, m+1:end), 'all') + sum(A(:, 1:m).*Z(:, 1:m), 'all');
term1 = 1/(p+1)*norm(Z(:, 1:m), 'fro')^2;
for ii = 1:p
   term1 = term1 + 2/(p+1-ii)*norm(Z(:, ii*m+1:(ii+1)*m), 'fro')^2;
end 
obj = 0.25*term1 + term2;
end