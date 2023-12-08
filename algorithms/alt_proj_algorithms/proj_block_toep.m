% Projects onto the set of block Toeplitz matrices.
function [out] = proj_block_toep(X, m, p)
    out = block_toep(mat2cell(block_averaging(X, m, p), m, m*ones(1, p+1))');
end