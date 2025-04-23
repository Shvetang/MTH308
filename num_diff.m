f = @(x, y) sin((x + y)^2);   % dy/dx = sin((x + y)^2)
x0 = 0;
y0 = -1;
n_values = [2, 6, 20, 63, 200, 632, 2000];

% Initialize table data storage
results = zeros(length(n_values), 3); % columns: n, rk4_error, euler_error

fprintf('%10s | %15s | %15s\n', 'n', 'RK4 Error', 'Euler Error');
fprintf('%s\n', repmat('-', 1, 45));

for k = 1:length(n_values)
    n = n_values(k);
    h = 4/n;

    [x_rk4, y_rk4] = rk4(f, x0, y0, h, n);
    rk4_diff = abs(y_rk4(end) - y_rk4(end-1));
    
    [x_euler, y_euler] = euler(f, x0, y0, h, n);
    euler_diff = abs(y_euler(end) - y_euler(end-1));
    
    results(k, :) = [n, rk4_diff, euler_diff];

    fprintf('%10d | %15.6e | %15.6e\n', n, rk4_diff, euler_diff);
end

% Compute for n = 200
n = 63;
h = 4/n;

[x_rk4, y_rk4] = rk4(f, x0, y0, h, n);
[x_euler, y_euler] = euler(f, x0, y0, h, n);


% Plotting the results
figure;
plot(x_rk4, y_rk4, 'b-', 'LineWidth', 1.5); hold on;
plot(x_euler, y_euler, 'r--', 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('RK4 vs Improved Euler Approximation (n = 200)');
legend('RK4', 'Improved Euler');
grid on;

% ---- Function Definitions ----

function [x, y] = rk4(f, x0, y0, h, n)
    x = zeros(1, n+1);
    y = zeros(1, n+1);
    x(1) = x0;
    y(1) = y0;

    for k = 1:n
        x_k = x(k);
        y_k = y(k);
        k1 = h * f(x_k, y_k);
        k2 = h * f(x_k + 0.5*h, y_k + 0.5*k1);
        k3 = h * f(x_k + 0.5*h, y_k + 0.5*k2);
        k4 = h * f(x_k + h, y_k + k3);
        x(k+1) = x_k + h;
        y(k+1) = y_k + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

function [x, y] = euler(f, x0, y0, h, n)
    x = zeros(1, n+1);
    y = zeros(1, n+1);
    x(1) = x0;
    y(1) = y0;

    for k = 1:n
        x_k = x(k);
        y_k = y(k);
        y_predict = y_k + h * f(x_k, y_k); % Euler predictor
        x(k+1) = x_k + h;
        y(k+1) = y_k + 0.5 * h * (f(x_k, y_k) + f(x(k+1), y_predict)); % Corrector
    end
end
