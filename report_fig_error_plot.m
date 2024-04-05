%% Geometric Multigrid
% In this file, we solve poisson equation -\delta u = f with dirichlet
% boundary condition. The domain we are using is [0,1] uniform grid for 1D 
% solver and [0,1] x [0,1] uniform grid for 2D solver.
clc;
close all;
%% setup
% Define mesh
% We require n to be an integer power of 2
n = 64; % thus delta x = delta y = h = 1/n. There will be (n-1) * (n-1) 
% grid points in the interior of the domain, and (n+1) * (n+1) grid points
% in the closure of the domain (i.e. include zeros on the boundary).
[X,Y] = meshgrid((1 : (n-1)) / n,(1 : (n-1)) / n);

% Define right-hand side f
%f = exp(-cos(4*X).^2 + exp(-sin(6*Y).^2)) - 3/2;
f = @(X,Y) (sqrt((X - 0.35).^2 + (Y - 0.6).^2) <= 0.1) * 1 + ...
           (sqrt((X - 0.8).^2 + (Y - 0.25).^2) <= 0.1) * 1;
b = reshape(f(X,Y),[],1);

% Define the discretized Laplacian operator
A = discrete_Laplacian2D(n);

% Thus the exact solution is
u_exact = A \ b;

%% Define parameters for solver 
% Define smoothing_method in {Gauss-Seidel, Jacobi, Relaxation}
smoothing_method = "Gauss-Seidel";
omega = 1;
max_iter = 30; % max iteration
nu1 = 2; % pre-smoothing step
nu2 = 2; % post-smoothing steps
gamma = 2; % V-cycle: gamma = 1, W-cycle: gamma = 2
tol = 1e-15;

%% solve & track error/convergence: W-cycle
iter = 1;
u_exact_matrix = reshape(u_exact, n-1, n-1);
error = zeros((n-1)^2, 1);
% iterations of multigrid cycle
%u_init = zeros((n-1)^2, 1);
u_init = 0.1 * rand((n-1)^2, 1);
u = u_init;
error(1) = norm(u_init - u_exact, Inf);
while iter <= max_iter && norm(u - u_exact, Inf) > tol
    fprintf('starting iteration %i...\n',iter);
    u = multigrid(@discrete_Laplacian2D, @restriction2D, b, n, gamma, u, nu1, nu2, smoothing_method, omega);
    error = reshape(u - u_exact, n-1, n-1);
    if rem(iter,2) == 1
        f1 = figure("Name", "Error After Cycle" + iter);
        surf(X, Y, error);
        title(sprintf("W-Cycle Error After Cycle %d",iter), 'FontSize', 20);
        exportgraphics(f1,"fig/w_cycle/MGW_Error_After_iter"+iter+".png","Resolution", 300, 'BackgroundColor','white');
    end
    iter = iter + 1;
end
fprintf('W cycle done. \n');
%% solve & track error/convergence: V-cycle
gamma = 1;
u = u_init;
error(1) = norm(u_init - u_exact, Inf);
iter = 1;
while iter <= max_iter && norm(u - u_exact, Inf) > tol
    fprintf('starting iteration %i...\n',iter);
    u = multigrid(@discrete_Laplacian2D, @restriction2D, b, n, gamma, u, nu1, nu2, smoothing_method, omega);
    u_init = u;
    error = reshape(u - u_exact, n-1, n-1);
    if rem(iter,2) == 1
        f1 = figure("Name", "Error After Cycle" + iter);
        surf(X, Y, error);
        title(sprintf("V-Cycle Error After Cycle %d",iter), 'FontSize', 20);
        exportgraphics(f1,"fig/v_cycle/MGV_Error_After_iter"+iter+".png","Resolution", 300, 'BackgroundColor','white');
    end
    iter = iter + 1;
end
fprintf('V cycle done. \n');