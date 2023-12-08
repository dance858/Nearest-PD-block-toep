function [X_feasible, X, obj, iter] = alt_proj_2(D, max_iter, m, p, tol, ref_value)
X = D;
P = zeros(m*(p+1), m*(p+1));
Q = zeros(m*(p+1), m*(p+1));
for k = 1:max_iter
   
    Y = proj_PSD_cone(X + P);
    P = X + P - Y;
    X = proj_block_toep(Y + Q, m, p);
    Q = Y + Q - X;
    
    if ~rem(k, 20)
        fprintf("Iteration: %i \n", k);
        X_feasible = X - min(0, min(eig(X)))*eye(m*(p+1)); 
        obj = norm(X_feasible - D, 'fro')^2;
        
        if (obj - ref_value)/ref_value <= tol
            break
        end
    end
end
iter = k;
end