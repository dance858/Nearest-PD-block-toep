function [output] = block_averaging(A, m, p)
% Averages the lower diagonal blocks of a matrix A.
output = zeros(m, m*(p+1));

%for block_row = 0:p
%    for block_col = 0:p-block_row
for block_col = 0:p
    for block_row = 0:p-block_col
        output(1:m, block_row*m+1:(block_row+1)*m) = ...
            output(1:m, block_row*m+1:(block_row+1)*m) + ...   
            A((block_row + block_col)*m + 1: (block_row + block_col + 1)*m, ...
                 block_col*m + 1: (block_col + 1)*m);
    end       %output(:, block_row*p + 1: block_row*(p+1)) = 1; 
end

for block = 0:p
   output(1:m, block*m + 1:(block+1)*m) = output(1:m, block*m + 1:(block+1)*m)/(p+1-block);
end

end