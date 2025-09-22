clear;
solver_params = struct();
solver_params.numerical_diff = 1;


[x_guess,exit_flag] = multi_newton(@object_dist, [1, 3], solver_params)

%run_simulation(theta,t_c)

function fval = object_dist(x) 
    % x = theta, t
    fval = projectile_traj(x(1),x(2)) - target_traj(x(2));
end