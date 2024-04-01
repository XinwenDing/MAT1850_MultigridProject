function u = multigrid(discrete_Laplacian, restriction, b, n, gamma, u_init, nu1, nu2, smoothing_method, omega)
    fprintf("n = %d\n multigrid starts, size A: ", n);
    A = discrete_Laplacian(n);
    disp(size(A));
    if n == 4
        u = A \ b;
    else
        R = restriction(n);
        P = 4 * R';
        % pre-smoothing
        fprintf("Pre-smoothing...\n");
        if smoothing_method == "Jacobi"
            u = Jacobi(A, b, u_init, nu1, omega);
        elseif smoothing_method == "Gauss-Seidel"
            u = GaussSeidel(A, b, u_init, nu1, omega);
        end

        % compute the defect
        res = b - A * u;

        % restrict the defect
        res = R * res;

        % compute an (approximate) solution v^tilde_{n-1} of the defect equation on Omega_{n/2}
        e_init = zeros((n/2 - 1)^2,1);	% start with zero
        for j = 1 : gamma
            fprintf("j = %d\n", j);
            e = multigrid(discrete_Laplacian, restriction, res, n/2, gamma, e_init, nu1, nu2, smoothing_method, omega);
        end

        % interpolate the correction
        e = P * e;

        % compute the corrected approximation on Omega_n
        u = u + e;

        % post-smoothing
        fprintf("Post-smoothing n = %d ...\n",n);
        if smoothing_method == "Jacobi"
            u = Jacobi(A, b, u, nu2, omega);
        elseif smoothing_method == "Gauss-Seidel"
            u = GaussSeidel(A, b, u, nu2, omega);
        end
    end
end