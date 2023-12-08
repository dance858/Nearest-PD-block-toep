function [s_prime, s_prime_prime] = s_diff(U, m, c, lambda, kernel)
% Evaluates the first and second derivative of the dual objective s.
% If kernel == 1, the function s is defined as
%       s(lambda) = log det(E'*[Toep(C^T) + lambda*I]^{-1} E) + c*lambda.
% If kernel == 2, the function s is defined as
%       s(lambda) = log det(E'*[Toep(C^T) + lambda*I]^{-1} E) + lambda^2/(ca).
% INPUT:
%       U: An upper triangular matrix U satisfying 
%         (Toep(C^T) + lambda*I)^{-1} = U*U';
%       m: dimension of the blocks

% Compute U*U'*E
temp1 = U*U(1:m, :)';

% Compute B = E'*U*U'*Toep(N^T)*U*U'E.
B = temp1'*temp1;

% Factor E'*(Toep(C^T) + lambda*I)^{-1} E = R'*R
R = chol(U(1:m, :)*U(1:m, :)');

first_deriv_matrix = (R\(R'\B));

% Compute first derivative.
if kernel == 1
    s_prime = c - trace(first_deriv_matrix);
else
    s_prime = 1/(2*c)*lambda - trace(first_deriv_matrix);
end


% Compute U'*Toep(N^T)*U*U'*E
temp3 = U'*temp1;

% Compute E'*U*U'*Toep(N^T)*U*U'*Toep(N^T)*U*U'*E.
temp4 = temp3'*temp3;

s_prime_prime = 2*trace((R\(R'\temp4))) - trace(first_deriv_matrix*first_deriv_matrix);

if kernel == 2
   s_prime_prime = s_prime_prime + 1/(2*c); 
end

end