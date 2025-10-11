function [global_trunc_error, num_evals] = calc_global_trunc_error(solver, rate_func_in, analytical_soln, tspan, h_ref)

    % Get the current state from the analytical solution
    X0 = analytical_soln(tspan(1)); % start from exact initial condition

    % Compute next numerical step
    switch solver
        case 'ForwardEuler'
           [t_list, X_list_num, ~, num_evals] = forward_euler_fixed_step(rate_func_in, tspan, X0, h_ref);
           X_list_ana = analytical_soln(t_list(end));
           global_trunc_error = norm(X_list_num(end, :)' - X_list_ana);

        case 'ExplicitMidpoint'
           [t_list, X_list_num, ~, num_evals] = explicit_midpoint_fixed_step(rate_func_in, tspan, X0, h_ref);
           X_list_ana = analytical_soln(t_list(end));
           global_trunc_error = norm(X_list_num(end, :)' - X_list_ana);

        otherwise
            warning('Invalid method');
            global_trunc_error = NaN;
            return;
    end
end
