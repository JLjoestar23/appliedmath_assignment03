% Compute the Jacobian matrix of a multivariate, vector-valued function 
% numerically using central finite differences.
%
% INPUTS:
% fun: a function handle that maps an input vector x (n×1) to an output 
% vector f(x) (m×1).
% x: nx1 column vector at which the Jacobian is to be approximated.
%
% OUTPUTS:
% J: approximated mxn Jacobian of fun evaluated at x. Each entry J(i,j) 
% approximates ∂f_i/∂x_j
% num_evals: the number of times the input function was called.
function [J, num_evals] = approximate_jacobian(fun, x)
    
    % initialize dimensions of J
    input_dim = length(x);
    output_dim = length(fun(x));

    J = zeros(output_dim, input_dim); % initialize empty J

    h = 1e-6; % finite difference step
    
    % using finite differences to evaluate each row of J
    for i=1:input_dim
        std_basis_vec = zeros(input_dim, 1);
        std_basis_vec(i) = h/2;
        J(:, i) = (fun(x + std_basis_vec) - fun(x - std_basis_vec)) / h;
    end
    
    num_evals = 2; % input function is called twice here
end