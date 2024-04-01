function [y, iter, err] = Jacobi_with_tol(A, b, y_init, max_iter, tol, true_soln)
    n = size(b,1);
    y = y_init;
    current_error = Inf;
    iter = 0;
    err = zeros(max_iter, 1);
    while iter < max_iter && current_error >= tol
        yprev = y;
        iter = iter + 1;
        for i = 1 : n
            idx_behind = 1 : (i-1);
            idx_ahead = (i+1) : n;
            y(i) = (b(i)- A(i,[idx_behind, idx_ahead]) * ...
                          yprev([idx_behind, idx_ahead])) / A(i,i);
        end
        % store error
        current_error = max(abs(true_soln - y));
        err(iter) = current_error;
    end
end
