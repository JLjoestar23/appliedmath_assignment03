%% Testing Forward Euler Fixed Step on test function 1

test_func_num = @rate_func01;
test_func_analytical = @solution01;

tspan = [0; 3*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = forward_euler_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('h = %.2f', h_ref(i)));
end

title('Accuracy of Forward Euler Integration for Varying Step Sizes');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Testing Forward Euler Fixed Step method on test function 2

test_func_num = @rate_func02;
test_func_analytical = @solution02;

tspan = [0; 3*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = forward_euler_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('h = %.2f', h_ref(i)));
end

title('Accuracy of Forward Euler Integration for Varying Step Sizes');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Explicit Midpoint fixed step test on test function 1

test_func_num = @rate_func01;
test_func_analytical = @solution01;

tspan = [0; 3*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = 1;
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
plot(t, test_func_analytical(t), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');

for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = explicit_midpoint_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list, 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('h = %.2f', h_ref(i)));
end

title('Accuracy of Explicit Midpoint Integration for Varying Step Sizes');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Testing Explicit Midpoint Method Fixed Step on test function 2

test_func_num = @rate_func02;
test_func_analytical = @solution02;

tspan = [0; 3*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = explicit_midpoint_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('h = %.2f', h_ref(i)));
end

title('Accuracy of Explicit Midpoint Integration for Varying Step Sizes');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Calculate and plot local truncation error for Forward Euler and Explicit Midpoint Method

method_list = {'ForwardEuler', 'ExplicitMidpoint'};
rate_func = @rate_func01;
analytical_soln = @solution01;
t_ref = 0.314;
num_trials = 1000;

local_trunc_error_analysis(method_list, rate_func, analytical_soln, t_ref, num_trials);

%% Calculate and plot global truncation error for Forward Euler and Explicit Midpoint Method

method_list = {'ForwardEuler', 'ExplicitMidpoint'};
rate_func = @rate_func01;
analytical_soln = @solution01;
tspan = [0, 2*pi];
num_trials = 200;

global_trunc_error_analysis(method_list, rate_func, analytical_soln, tspan, num_trials);

%% Testing Implicit Euler Fixed Step method on test function 01

test_func_num = @rate_func01;
test_func_analytical = @solution01;
step_func = @implicit_euler_step;

tspan = [0; 3*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = fixed_step_integration(test_func_num, step_func, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('h = %.2f', h_ref(i)));
end

title('Accuracy of Implicit Euler Integration for Varying Step Sizes');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Testing Implicit Euler Fixed Step method on test function 02

test_func_num = @rate_func02;
test_func_analytical = @solution02;
step_func = @implicit_euler_step;

tspan = [0; 3*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = fixed_step_integration(test_func_num, step_func, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('h = %.2f', h_ref(i)));
end

title('Accuracy of Implicit Euler Integration for Varying Step Sizes');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Testing Implicit Midpoint Fixed Step method on test function 02

test_func_num = @rate_func02;
test_func_analytical = @solution02;
step_func = @implicit_midpoint_step;

tspan = [0; 3*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = fixed_step_integration(test_func_num, step_func, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', sprintf('h = %.2f', h_ref(i)));
end

title('Accuracy of Implicit Euler Integration for Varying Step Sizes');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Local truncation error comparison bewteen explicit and implicit methods

method_list = {'ForwardEuler', 'ExplicitMidpoint', 'ImplicitEuler', 'ImplicitMidpoint'};
rate_func = @rate_func01;
analytical_soln = @solution01;
t_ref = 0.314;
num_trials = 1000;

local_trunc_error_analysis_modified(method_list, rate_func, analytical_soln, t_ref, num_trials);

%% Global truncation error comparison bewteen explicit and implicit methods

method_list = {'ForwardEuler', 'ExplicitMidpoint', 'ImplicitEuler', 'ImplicitMidpoint'};
rate_func = @rate_func01;
analytical_soln = @solution01;
tspan = [0, 2*pi];
num_trials = 500;

global_trunc_error_analysis_modified(method_list, rate_func, analytical_soln, tspan, num_trials);