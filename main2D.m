%% Geometric Multigrid
% In this file, we solve poisson equation -\delta u = f with dirichlet
% boundary condition. The domain we are using is [0,1] uniform grid for 1D 
% solver and [0,1] x [0,1] uniform grid for 2D solver.
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
max_iter = 50; % max iteration
nu1 = 2; % pre-smoothing step
nu2 = 2; % post-smoothing steps
gamma = 2; % V-cycle: gamma = 1, W-cycle: gamma = 2
tol = 1e-12;

%% solve & track error/convergence: W-cycle
iter = 1;
u_exact_matrix = reshape(u_exact, n-1, n-1);
error = zeros((n-1)^2, 1);
% iterations of multigrid cycle
%u_init = 0.1 * rand((n-1)^2, 1);
u_init = zeros((n-1)^2, 1);
u = u_init;
error(1) = norm(u_init - u_exact, Inf);
while iter <= max_iter && norm(u - u_exact, Inf) > tol
    u = multigrid(@discrete_Laplacian2D, @restriction2D, b, n, gamma, u, nu1, nu2, smoothing_method, omega);
    error = reshape(u - u_exact, n-1, n-1);
    iter = iter + 1;
end

%% Enforce BC and plot solution
% define meshgrid with boundary data point
[X_BC,Y_BC] = meshgrid(0:1/n:1, 0:1/n:1);
% plot solution from multigrid solver
u = padarray(reshape(u, n-1, n-1), [1, 1], 0);
figure("Name", "Multigrid Solution 2D")
surf(X_BC, Y_BC, u);
title("Multigrid Solution of Heat Equation");
% plot exact solution 
u_exact = padarray(reshape(u_exact, n-1, n-1), [1, 1], 0);
figure("Name", "Exact Solution 2D")
surf(X_BC, Y_BC, u_exact);
title("Exact Solution of Heat Equation");
