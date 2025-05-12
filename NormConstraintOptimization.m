clear; clc;
%% Set Up the Initial Guess
% Generate initial guesses for x and y
N = 3; % Example value for N
x0 = randn(N-1, 1);
y0 = randn(N-1, 1);

% Normalize x0 and y0 to satisfy the norm constraints
x0 = x0 / norm(x0);
y0 = y0 / norm(y0);

% Combine x0 and y0 into a single vector
vars0 = [x0; y0];

%% Call the Optimization Function
% Optimization options
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');


% Define anonymous functions to pass additional parameters
objFun = @(vars) objectiveFunction(vars, N);
conFun = @(vars) constraintFunction(vars, N);

% Run the optimization
[vars_min, C_min] = fmincon(objFun, vars0, [], [], [], [], [], [], conFun, options);

% Extract optimized x and y
x_min = vars_min(1:N-1);
y_min = vars_min(N:2*N-2);

% Display the results
fprintf('Minimum value of C: %f\n', C_min);

%% Verify the solution
% Check norms
norm_x = norm(x_min);
norm_y = norm(y_min);
fprintf('Norm of x_min: %f (should be 1)\n', norm_x);
fprintf('Norm of y_min: %f (should be 1)\n', norm_y);

% Check inner product constraint
rho_min = dot(x_min, y_min);
fprintf('Inner product squared: %f (should be less than 1)\n', rho_min^2);

%% Plot

figure;
p = fplot(@(t) 1*sin(t), @(t) 1*cos(t), [0, 2*pi], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5);
hold on;
plot([-1, 1], [0, 0], 'k', 'LineStyle', ':', 'LineWidth', 0.5); % x축 기준선
plot([0, 0], [-1, 1], 'k', 'LineStyle', ':', 'LineWidth', 0.5); % y축 기준선
scatter(x_min, y_min, 'b', 'filled');
xlabel('x');
ylabel('y');
title('Optimal anchor placement');
grid on; % 그리드를 추가하여 보기 쉽게 합니다
scatter(0,0, 'b','filled');
hold off;

%%
function [c, ceq] = constraintFunction(vars, N)
    % Extract x and y from vars
    x = vars(1:N-1);
    y = vars(N:2*N-2);
    
    % Compute norms
    norm_x_sq = sum(x.^2);
    norm_y_sq = sum(y.^2);
    
    % Equality constraints 
    % norm_x_sq = norm_y_sq = 1, 
    ceq = [norm_x_sq - 1;
           norm_y_sq - 1];
    % Inequality constraint: 
    c = [];
end
%%
function C = objectiveFunction(vars, N)
    % Extract x and y from vars(= vars consists of x and y)
    x = vars(1:N-1);
    y = vars(N:2*N-2);
    
    % Compute sums
    sum_x = sum(x);
    sum_y = sum(y);
    norm_x = norm(x);
    norm_y = norm(y);
    
    % Compute inner product
    rho = dot(x, y)/(norm_x*norm_y);
    
    % Compute numerator and denominator
    numerator = (sum_x - rho * sum_y)^2 + (sum_y - rho * sum_x)^2 + 2 * (1 - rho^2);
    denominator = (1 - rho^2)^2;
    % Compute C
    C = numerator / denominator;
end
