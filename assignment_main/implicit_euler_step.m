% This function computes the value of X at the next time step using the 
% Implicit Euler approximation
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
function [XB, num_evals] = implicit_euler_step(rate_func_in, t, XA, h)
    % initializing solver params
    solver_params.approx_j = 1;
    
    % setting up implicit equation to numerically solve for XB
    solve_for_XB_implicit = @(XB) XA + h*rate_func_in(t + h, XB) - XB;
    
    % numerically solving for XB
    [XB, num_evals] = multivariate_newton_solver(solve_for_XB_implicit, XA, solver_params);

end