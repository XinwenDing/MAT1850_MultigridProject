function R = restriction2D(n)
    R = restriction1D(n);
    R = kron(R,R);
end
