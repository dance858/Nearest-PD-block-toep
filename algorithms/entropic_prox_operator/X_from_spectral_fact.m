function [X] = X_from_spectral_fact(B, m, p)
% X = D(B*B')'. But weird convention on B.

X = zeros(m, m*(p+1));

for k = 0:p
   for i = 0:p-k
        X(:, k*m+1:(k+1)*m) = X(:, k*m+1:(k+1)*m) + B(i*m+1:(i+1)*m, :)*B((i+k)*m+1:(i+k+1)*m, :)';
   end
end

end