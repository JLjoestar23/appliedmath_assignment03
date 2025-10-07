%% Forward Euler fixed step test
test_func_num = @rate_func01;
test_func_analytical = @solution01;

tspan = [0; 4*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = 1;
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
plot(t, test_func_analytical(t), 'b-', 'LineWidth', 1.5);

for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = forward_euler_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list, '.', 'MarkerSize', 10);
end

title('Comparing Forward Euler and Analytical Solutions for Test Function 01');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;

%% Explicit Midpoint fixed step test
test_func_num = @rate_func01;
test_func_analytical = @solution01;

tspan = [0; 4*pi];
t = linspace(tspan(1), tspan(2), 500);
X0 = 1;
h_ref = (0.1: 0.1: 0.4);

figure();
hold on;
plot(t, test_func_analytical(t), 'b-', 'LineWidth', 1.5);

for i=1:length(h_ref)
    [t_list, X_list, h_avg_1, num_evals_1] = explicit_midpoint_fixed_step(test_func_num, tspan, X0, h_ref(i));
    plot(t_list, X_list, '.', 'MarkerSize', 10);
end

title('Comparing Explicit Midpoint and Analytical Solutions for Test Function 01');
legend();
grid on;
xlim([t_list(1) t_list(end)])
hold off;