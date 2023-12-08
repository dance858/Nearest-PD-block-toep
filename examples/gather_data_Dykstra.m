clear; clc;
addpath('../algorithms/alt_proj_algorithms')
addpath('../algorithms/entropic_prox_operator')

%% Settings for experiment. Note that 
m = 20;
all_p = [20, 40, 60, 80];
num_runs = 5;

%% Settings for Dykstras 
max_iter = 4000; tol = 1e-3;

for p = all_p
    fprintf("Current p: %i\n", p);
    for run = 2:num_runs
        % Read in bound from IGA to get a "fair" termination criteria
        % Make sure to run gather_data_IGA before this.
        filename = "data_IGA/m=" + num2str(m) + "/BT_m=" + num2str(m) + "_p=" + ...
                    num2str(p) + "_run=" + num2str(run) + "_c=0.001_L=5.mat";
        load(filename, "BT", 'UB_primal');
        
        start = tic;
        [X_feasible, X, obj, iter]  = alt_proj_2(BT, max_iter, m, p, tol, UB_primal);
        solve_time = toc(start);
        
        filename = "data_dykstra/m=" + num2str(m) + "/BT_m=" + num2str(m) + "_p=" ...
                  + num2str(p) + "_tol=" + num2str(tol) + "_run="  ...
                  + num2str(run) + ".mat";
        
        save(filename);
    end
end


