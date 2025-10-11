function [local_trunc_error, local_h_diff] = calc_local_trunc_error(solver, rate_func_in, analytical_soln, t, h)
    
    % Get the current state from the analytical solution
    XA_num = analytical_soln(t); % start from exact initial condition

    % Analytical next step
    XB_ana = analytical_soln(t + h);
    switch solver
        case 'ForwardEuler'
            % Compute derivative at current state
            dXdt = rate_func_in(t, XA_num);

            % Numerical next step
            XB_num = XA_num + h*dXdt;

        case 'ExplicitMidpoint'
            % Derivative at current state
            dXdt_1 = rate_func_in(t, XA_num);

            % Midpoint estimate
            XB_half = XA_num + (h/2) * dXdt_1;

            % Derivative at midpoint
            dXdt_2 = rate_func_in(t + h/2, XB_half);

            % Numerical next step
            XB_num = XA_num + h*dXdt_2;

        otherwise
            warning('Invalid method');
            local_trunc_error = NaN;
            return;
    end

    % calculating local truncation error
    local_trunc_error = norm(XB_num - XB_ana);

    % calcaulting local timestep difference based on analytical sol'n
    local_h_diff = norm(analytical_soln(t+h) - analytical_soln(t));
end
