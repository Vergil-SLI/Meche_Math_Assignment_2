theta = pi/4;

solver_params = struct();
solver_params.numerical_diff = 1;
object_dist = @(t) (projectile_traj(theta,t) - target_traj(t));
[x_guess,exit_flag] = multi_newton(object_dist, 8, solver_params);

%run_simulation(theta,t_c)