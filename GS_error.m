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