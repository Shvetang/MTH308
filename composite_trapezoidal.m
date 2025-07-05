function I = compositetrapezoidal(f, a, b, n)
    % composite_trapezoidal - Composite Trapezoidal Rule
    % Inputs:
    %   f - function handle of the integrand
    %   a - lower limit of integration
    %   b - upper limit of integration
    %   n - number of subintervals
    % Output:
    %   I - approximate value of the integral

    h = (b - a) / n;                    % Step size
    x = a:h:b;                          % n+1 points
    y = f(x);                           % Function values at those points

    I = h * (0.5*y(1) + sum(y(2:end-1)) + 0.5*y(end));  % Trapezoidal rule
end

f = @(x) sin(x);              % Function to integrate
a = 0;
b = pi;
n = 100;                      % Number of subintervals

I = compositetrapezoidal(f, a, b, n);
fprintf("The Composite Trapezoidal Method is a numerical integration technique used to approximate the definite integral of a function over an interval [a,b]. \nIt works by dividing the interval into smaller subintervals, applying the Trapezoidal Rule (where the function is approximated by the average of its values at both endpoints) \non each, and summing up the areas of the trapezoids.")
fprintf('\nApproximate integral: %.6f\n', I);
