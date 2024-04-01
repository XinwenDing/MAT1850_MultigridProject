function y = Jacobi(A, b, y_init, max_iter, omega)
    n = size(b,1);
    y = y_init;
    iter = 0;
    while iter < max_iter
        yprev = y;
        iter = iter + 1;
        for i = 1 : n
            idx_behind = 1 : (i-1);
            idx_ahead = (i+1) : n;
            y(i) = (b(i)- A(i,[idx_behind, idx_ahead]) * ...
                          yprev([idx_behind, idx_ahead])) / A(i,i);
        end
        if omega ~= 1
            % Note: if omega == 1, then y = 0 * yprev + y
            y = (1 - omega) * yprev + omega * y;
        end
    end
end
