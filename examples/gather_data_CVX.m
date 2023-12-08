clear; clc;
m = 5;
all_p = [20];
num_runs = 2;
%solver = "scs";
solver = "SDPT3";

for p = all_p
    for run = 1:num_runs
        filename = "data/m=" + num2str(m) + "/BT_m=" + num2str(m) + "_p=" + ...
                    num2str(p) + "_run=" + num2str(run) + ".mat";
        load(filename, "BT");
                
        A_CVX = BT(:, 1:m);         % first block column
        block_toep_A = block_toep(mat2cell(A_CVX, m*ones(1, p+1), m)');
        [TX_primal, cvx_optval_primal, solve_time] = solve_primal_with_cvx(A_CVX, m, p, solver);
        
        fprintf("solve_time: %f \n", solve_time);
        filename = "data_cvx/m=" + num2str(m) + "/BT_m=" + num2str(m) + "_p=" + num2str(p) + "_run=" + ...
                    num2str(run) + "_solver=" + solver + ".mat";
        save(filename); 
    end
end


