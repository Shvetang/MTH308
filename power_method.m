function [lambda, x, iter, history] = power_method(A, x0, max_iter)
    % POWER_METHOD finds the dominant eigenvalue and eigenvector using the power method
    % Inputs:
    %   A         - square matrix
    %   x0        - initial guess vector
    %   max_iter  - maximum number of iterations (optional)
    % Outputs:
    %   lambda    - estimated dominant eigenvalue
    %   x         - corresponding eigenvector
    %   iter      - number of iterations

    if nargin < 3
        max_iter = 10; % default max iterations
    end

    x = x0 / norm(x0, inf); % normalize initial guess
    iter = 0;

    % initialize storage for history
   n = length(x0);

   varNames = ["iter_no", compose("y%d", 1:n), "mu", compose("x%d", 1:n)];
   varNames = string(varNames); % Ensure all are strings

   history = table('Size', [0, 2*n + 2], ...
                'VariableTypes', repmat("double", 1, 2*n + 2), ...
                'VariableNames', varNames);

    while iter < max_iter
        iter = iter + 1;

        y = A * x;
        [~, i] = max(abs(y));   % find index of maximum absolute value
        mu = y(i);              % corresponding value for scaling
        x_new = y / mu;

        % Save to history table
        history(iter, :) = array2table([iter, y(:)', mu, x_new(:)']);

        x = x_new;
    end

    lambda = mu;
    x = x_new; % final normalized eigenvector
end

A = [2 2 2; 2/3 5/3 5/3; 1 5/2 11/2];
x0 = [1; 0; 0];

[lambda, x, iter, history] = power_method(A, x0);

disp("Iteration history:")
disp(history)

fprintf('Dominant eigenvalue: %.6f\n', lambda);
fprintf('Corresponding eigenvector:');
disp(x');
