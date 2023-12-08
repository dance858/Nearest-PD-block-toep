function [grad] = grad_f_conj(Z, A, m, p)
grad = A;
for k = 0:p 
       grad(:, k*m + 1:(k+1)*m) = grad(:, k*m + 1:(k+1)*m) ...
                                + 0.5/(p+1-k)*Z(:, k*m+1:(k+1)*m);
end
end