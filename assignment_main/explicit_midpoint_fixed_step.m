% Runs numerical integration using explicit midpoint approximation
% 
% INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
% tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
% X0: the vector describing the initial conditions, X(t_start)
% h_ref: the desired value of the average step size (not the actual value)
% 
% OUTPUTS:
% t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
% X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
% h_avg: the average step size
% num_evals: total number of calls made to rate_func_in during the integration
function [t_list, X_list, h_avg, num_evals] = explicit_midpoint_fixed_step(rate_func_in, tspan, X0, h_ref)
    % Calculate number of steps so that step size h is as close to h_ref as possible
    N = floor((tspan(2) - tspan(1)) / h_ref);

    % Actual step size
    h_avg = (tspan(2) - tspan(1)) / N;

    % Time vector
    t_list = linspace(tspan(1), tspan(2), N+1);
    
    % Preallocate solution array (each row = one time step)
    X_list = zeros(N+1, length(X0));
    X_list(1,:) = X0;

    num_evals = 0; % initialize function call counter

    for i = 1:N
        % Evaluate next step using explicit midpoint
        [XB, evals] = explicit_midpoint_step(rate_func_in, t_list(i), X_list(i,:)', h_avg);
        
        % Store result (convert to row vector for X_list)
        X_list(i+1,:) = XB(:)';

        % Accumulate total number of evaluations
        num_evals = num_evals + evals;
    end
end

% This function computes the value of X at the next time step using the 
% explicit midpoint approximation
%
% INPUTS:
% rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
% t: the value of time at the current step
% XA: the value of X(t)
% h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%
% OUTPUTS:
% XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
% num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB, num_evals] = explicit_midpoint_step(rate_func_in, t, XA, h)
    % Step 1: derivative at the start
    dXdt_1 = rate_func_in(t, XA);

    % Step 2: estimate midpoint value
    XB_half = XA + (h/2) * dXdt_1;

    % Step 3: derivative at midpoint
    dXdt_2 = rate_func_in(t + h/2, XB_half);

    % Step 4: use midpoint derivative for final update
    XB = XA + h * dXdt_2;

    % Explicit midpoint requires 2 function calls per step
    num_evals = 2;
end