clear; close all; clc;
%% Set Up the Initial Guess
% Generate initial guesses for x and y
N = 4; % Number of anchors (excluding reference)
x0 = randn(N, 1);
y0 = randn(N, 1);

% Combine x0 and y0 into a single decision vector
vars0 = [x0; y0];

%% Call the Optimization Function
% Choose SQP algorithm for fmincon
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');

% Precompute Delaunay edges on the initial guess (fixed neighbor pairs)
DT0 = delaunayTriangulation(x0, y0);
E   = edges(DT0);             % M×2 array of neighbor index pairs

% Nonlinear constraint function handle uses fixed edge list E
conFun = @(vars) constraintFunction(vars, N, E);

% Objective function handle
objFun = @(vars) objectiveFunction(vars, N);

%% Boundary Constraints for First Anchor
% No bounds except to fix the first anchor at (1,0)
lb = -inf(2*N, 1);
ub =  inf(2*N, 1);
lb(1)     = 1;  % Fix x(1) = 1
ub(1)     = 1;
lb(N+1)   = 0;  % Fix y(1) = 0
ub(N+1)   = 0;

% Run constrained optimization
[vars_min, C_min] = fmincon(objFun, vars0, [], [], [], [], lb, ub, conFun, options);

% Extract optimized x and y
x_min = vars_min(1:N);
y_min = vars_min(N+1:2*N);

%% Display Results
fprintf('Dot product (x, y): %f\n', dot(x_min, y_min));
fprintf('Sum(x), Sum(y): %f, %f\n', sum(x_min), sum(y_min));
fprintf('Max radius: %f\n', max(sqrt(x_min.^2 + y_min.^2)));
fprintf('Minimum cost C: %f\n', C_min);
disp('Optimized coordinates [x, y]:');
disp([x_min, y_min]);

%% Plot Results
figure;
% Plot unit circle
fplot(@(t) cos(t), @(t) sin(t), [0, 2*pi], 'k:', 'LineWidth', 0.5);
axis equal; hold on;
% Scatter optimized points
scatter(x_min, y_min, 'b', 'filled');
% Draw axes
plot([-1,1],[0,0],'k:','LineWidth',0.5);
plot([0,0],[-1,1],'k:','LineWidth',0.5);
xlim([-1.1, 1.1]); ylim([-1.1, 1.1]);
xlabel('x'); ylabel('y');
title('Optimal Anchor Placement on Unit Circle');
hold off;

% Display angles for each anchor
disp('Angles theta_n:');
disp(atan2(y_min, x_min));



%% Objective Function
function C = objectiveFunction(vars, N)
    % Split vars into x and y
    x = vars(1:N);
    y = vars(N+1:2*N);
    
    % Build difference matrix A for anchors 2..N relative to anchor 1
    anchor_matrix = [x, y];
    A = zeros(N-1, 2);
    for i = 2:N
        A(i-1, :) = anchor_matrix(i, :) - anchor_matrix(1, :);
    end
    
    % Compute beamforming cost components
    barx = A(:,1);
    bary = A(:,2);
    sx   = sum(barx);
    sy   = sum(bary);
    p    = sum(barx .* bary);
    nx2  = sum(barx.^2);
    ny2  = sum(bary.^2);
    
    Vx = (sx*ny2 - sy*p)^2;
    Vy = (sy*nx2 - sx*p)^2;
    for i = 1:N-1
        Vx = Vx + (barx(i)*ny2 - bary(i)*p)^2;
        Vy = Vy + (bary(i)*nx2 - barx(i)*p)^2;
    end
    denom = (nx2*ny2 - p^2)^2;
    Vx = Vx / denom;
    Vy = Vy / denom;
    
    % Total cost
    C = Vx + Vy;
end

%% Constraint Function
function [c, ceq] = constraintFunction(vars, N, E)
    % Split vars into x and y
    x = vars(1:N);
    y = vars(N+1:2*N);
    
    % Equality constraint: points lie on unit circle
    ceq = x.^2 + y.^2 - 1;
    
    % Inequality constraints: maintain minimum distance along Delaunay edges
    eps = 1e-1;
    M   = size(E,1);
    c   = zeros(M,1);
    for k = 1:M
        i = E(k,1);
        j = E(k,2);
        distSq = (x(j)-x(i))^2 + (y(j)-y(i))^2;
        c(k) = eps^2 - distSq;
    end
end
