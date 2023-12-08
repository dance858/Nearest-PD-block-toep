% Projects onto the set of positive semidefinite matrices.
function [out] = proj_PSD_cone(X)
   [u, v] = eig(X);
   out = u*max(0, v)*u';
end