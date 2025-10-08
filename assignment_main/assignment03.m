%% Testing Forward Euler Fixed Step on test function 1
test_func_num = @rate_func02;
test_func_analytical = @solution02;

tspan = [0; 2*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2);
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = forward_euler_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1);
end

title('Comparing Forward Euler and Analytical Solutions for Test Function 01');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Testing Forward Euler Fixed Step method on test function 2
test_func_num = @rate_func02;
test_func_analytical = @solution02;

tspan = [0; 2*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2);
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = forward_euler_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1);
end

title('Comparing Forward Euler and Analytical Solutions for Test Function 02');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Explicit Midpoint fixed step test
test_func_num = @rate_func01;
test_func_analytical = @solution01;

tspan = [0; 2*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = 1;
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
plot(t, test_func_analytical(t), 'b-', 'LineWidth', 2);

for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = explicit_midpoint_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list, 'o-', 'MarkerSize', 4, 'LineWidth', 1);
end

title('Comparing Explicit Midpoint and Analytical Solutions for Test Function 01');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Testing Explicit Midpoint Method Fixed Step on test function 2
test_func_num = @rate_func02;
test_func_analytical = @solution02;

tspan = [0; 2*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = [1, 0];
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
val = test_func_analytical(t);
plot(t, val(1, :), 'b-', 'LineWidth', 2);
for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = explicit_midpoint_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list(:, 1), 'o-', 'MarkerSize', 4, 'LineWidth', 1);
end

title('Comparing Forward Euler and Analytical Solutions for Test Function 02');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Calculate local truncation error

% initializing function handles
rate_func = @rate_func01;
analytical_soln = @solution01;

% arbitrary time to analyze truncation error at
t_ref = 0.314;

% number of error trials
num_trials = 1000;

% initialize list to plot errors
fe_local_truncation_errors = zeros(1, num_trials);
mp_local_truncation_errors = zeros(1, num_trials);

% initialize list of error values
%h_list = linspace(10e-6, 10, num_trials);
h_list = logspace(-6, 1, num_trials);

% sampling local truncation errors from Forward Euler and Midpoint
for i=1:num_trials
    fe_local_truncation_errors(i) = calc_local_trunc_error('ForwardEuler', rate_func, analytical_soln, t_ref, h_list(i));
    mp_local_truncation_errors(i) = calc_local_trunc_error('MidpointMethod', rate_func, analytical_soln, t_ref, h_list(i));
end

% line fit in log scale
start_i = 200;
end_i = 800;
logh = log10(h_list(start_i:end_i));
logerr_fe = log10(fe_local_truncation_errors(start_i:end_i));
logerr_mp = log10(mp_local_truncation_errors(start_i:end_i));

% Perform linear regression (polyfit on log–log data)
p_fe = polyfit(logh, logerr_fe, 1);
p_mp = polyfit(logh, logerr_mp, 1);

% Get slopes (order of accuracy)
slope_fe = p_fe(1);
slope_mp = p_mp(1);

% Generate fit lines for plotting
fit_fe = polyval(p_fe, logh);
fit_mp = polyval(p_mp, logh);

% plot the data
loglog(h_list, fe_local_truncation_errors, 'b.', 'MarkerSize', 5, 'DisplayName', 'Forward Euler');
hold on;
loglog(h_list, mp_local_truncation_errors, 'r.', 'MarkerSize', 5, 'DisplayName', 'Midpoint Method');
loglog(10.^logh, 10.^fit_fe, 'y--', 'LineWidth', 2, 'DisplayName', sprintf('FE Fit (slope = %.2f)', slope_fe));
loglog(10.^logh, 10.^fit_mp, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('MP Fit (slope = %.2f)', slope_mp));
title('Local Truncation Error vs. Step Size');
xlabel('Timestep (s)');
ylabel('Local Truncation Error');
legend('Location','northwest');
axis square;
grid on;
hold off;

%% Calculate global truncation error

% initializing function handles
rate_func = @rate_func01;
analytical_soln = @solution01;

% arbitrary time to analyze truncation error at
tspan = [0, 2*pi];

% number of error trials
num_trials = 250;

% initialize list to plot errors
fe_global_truncation_errors = zeros(1, num_trials);
mp_global_truncation_errors = zeros(1, num_trials);

fe_num_evals = zeros(1, num_trials);
mp_num_evals = zeros(1, num_trials);

% initialize list of error values
%h_list = linspace(10e-6, 10, num_trials);
h_list = logspace(-6, 1, num_trials);

% sampling local truncation errors from Forward Euler and Midpoint
for i=1:num_trials
    [fe_global_truncation_errors(i), fe_num_evals(i)] = calc_global_trunc_error('ForwardEuler', rate_func, analytical_soln, tspan, h_list(i));
    [mp_global_truncation_errors(i), mp_num_evals(i)] = calc_global_trunc_error('MidpointMethod', rate_func, analytical_soln, tspan, h_list(i));
end

% Generating plot for global truncation error vs. step size
% line fit in log scale
start_i = 50;
end_i = 150;
logh = log10(h_list(start_i:end_i));
logerr_fe = log10(fe_global_truncation_errors(start_i:end_i));
logerr_mp = log10(mp_global_truncation_errors(start_i:end_i));

% perform linear regression (polyfit on log–log data)
p_fe_1 = polyfit(logh, logerr_fe, 1);
p_mp_1 = polyfit(logh, logerr_mp, 1);

% get slopes (order of accuracy)
slope_fe_1 = p_fe_1(1);
slope_mp_1 = p_mp_1(1);

% generate fit lines for plotting
fit_fe_1 = polyval(p_fe_1, logh);
fit_mp_1 = polyval(p_mp_1, logh);

% plot the data
figure();
loglog(h_list, fe_global_truncation_errors, 'b.', 'MarkerSize', 5, 'DisplayName', 'Forward Euler');
hold on;
loglog(h_list, mp_global_truncation_errors, 'r.', 'MarkerSize', 5, 'DisplayName', 'Midpoint Method');
loglog(10.^logh, 10.^fit_fe_1, 'y--', 'LineWidth', 2, 'DisplayName', sprintf('FE Fit (slope = %.2f)', slope_fe_1));
loglog(10.^logh, 10.^fit_mp_1, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('MP Fit (slope = %.2f)', slope_mp_1));
title('Global Truncation Error vs. Step Size');
xlabel('Timestep (s)');
ylabel('Global Truncation Error');
legend('Location','northwest');
axis square;
grid on;
hold off;

% generating plot for global truncation error vs. number of function calls
% line fit in log scale
log_num_evals_1 = log10(fe_num_evals(start_i:end_i));
log_num_evals_2 = log10(mp_num_evals(start_i:end_i));

% perform linear regression (polyfit on log–log data)
p_fe_2 = polyfit(log_num_evals_1, logerr_fe, 1);
p_mp_2 = polyfit(log_num_evals_2, logerr_mp, 1);

% get slopes (order of accuracy)
slope_fe_2 = p_fe_2(1);
slope_mp_2 = p_mp_2(1);

% generate fit lines for plotting
fit_fe_2 = polyval(p_fe_2, log_num_evals_1);
fit_mp_2 = polyval(p_mp_2, log_num_evals_2);

% plot the data
figure();
loglog(fe_num_evals, fe_global_truncation_errors, 'b.', 'MarkerSize', 5, 'DisplayName', 'Forward Euler');
hold on;
loglog(mp_num_evals, mp_global_truncation_errors, 'r.', 'MarkerSize', 5, 'DisplayName', 'Midpoint Method');
loglog(10.^log_num_evals_1, 10.^fit_fe_2, 'y--', 'LineWidth', 2, 'DisplayName', sprintf('FE Fit (slope = %.2f)', slope_fe_2));
loglog(10.^log_num_evals_2, 10.^fit_mp_2, 'g--', 'LineWidth', 2, 'DisplayName', sprintf('MP Fit (slope = %.2f)', slope_mp_2));
title('Global Truncation Error vs. Number of Function Calls');
xlabel('Number of Function Calls');
ylabel('Global Truncation Error');
legend('Location','northeast');
axis square;
grid on;
hold off;

%% 