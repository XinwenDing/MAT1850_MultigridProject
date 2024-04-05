clc;
close all;
%% setup
% Define mesh
% We require n to be an integer power of 2
n = 32; % thus delta x = delta y = h = 1/n. There will be (n-1) * (n-1) 
% grid points in the interior of the domain, and (n+1) * (n+1) grid points
% in the closure of the domain (i.e. include zeros on the boundary).
[X,Y] = meshgrid((1 : (n-1)) / n,(1 : (n-1)) / n);

% Define right-hand side f
%f = exp(-cos(4*X).^2 + exp(-sin(6*Y).^2)) - 3/2;
f = @(X,Y) (sqrt((X - 0.35).^2 + (Y - 0.6).^2) <= 0.1) * 1 + ...
           (sqrt((X - 0.8).^2 + (Y - 0.25).^2) <= 0.1) * 1;
f = @(X,Y) exp(-cos(4*X).^2 + exp(-sin(6*Y).^2)) - 1.6;
b = reshape(f(X,Y),[],1);

% Define the discretized Laplacian operator
A = discrete_Laplacian2D(n);

% Thus the exact solution is
u_exact = A \ b;

%% Gauss-Seidel Error
tol = 1e-12;
%u_init = zeros((n-1)^2, 1);
u_init = 0.1 * rand((n-1)^2, 1);
max_iter1 = 1;
[u1, ~, err1] = GaussSeidel_with_tol(A, b, u_init, max_iter1, tol, u_exact);
max_iter2 = 2;
[u2, ~, err2] = GaussSeidel_with_tol(A, b, u_init, max_iter2, tol, u_exact);
max_iter3 = 3;
[u3, ~, err3] = GaussSeidel_with_tol(A, b, u_init, max_iter3, tol, u_exact);
u1 = reshape(u1, n-1, n-1);
u2 = reshape(u2, n-1, n-1);
u3 = reshape(u3, n-1, n-1);
u_exact = reshape(u_exact, n-1, n-1);
u_init = reshape(u_init, n-1, n-1);

%% plot error
%[X_BC,Y_BC] = meshgrid(0:1/n:1, 0:1/n:1);
[X_BC,Y_BC] = meshgrid((1 : (n-1)) / n,(1 : (n-1)) / n);
f0 = figure("Name", "Init error");
surf(X_BC, Y_BC, u_init - u_exact);
zlim([0 0.15])
title("Initial Error of GS", 'FontSize', 20);
f1 = figure("Name", "u1 error");
surf(X_BC, Y_BC, u1 - u_exact);
zlim([0 0.15])
title(sprintf("GS Error After %d Iteration", max_iter1), 'FontSize', 20);
f2 = figure("Name", "u2 error");
surf(X_BC, Y_BC, u2 - u_exact);
zlim([0 0.15])
title(sprintf("GS Error After %d Iterations", max_iter2), 'FontSize', 20);
f3 = figure("Name", "u3 error");
surf(X_BC, Y_BC, u3 - u_exact);
zlim([0 0.15])
title(sprintf("GS Error After %d Iterations", max_iter3), 'FontSize', 20);
% f4 = figure("Name", "u exact");
% surf(X_BC, Y_BC, u_exact);
% title(sprintf("exact solution"));

exportgraphics(f0,'fig/Initial_Error_of_GS.png',"Resolution", 300, 'BackgroundColor','white');
exportgraphics(f1,'fig/GS_Error_After_iter1.png',"Resolution", 300, 'BackgroundColor','white');
exportgraphics(f2,'fig/GS_Error_After_iter2.png',"Resolution", 300, 'BackgroundColor','white');
exportgraphics(f3,'fig/GS_Error_After_iter3.png',"Resolution", 300, 'BackgroundColor','white');

