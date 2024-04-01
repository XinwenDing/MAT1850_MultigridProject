function A = discrete_Laplacian2D(n)
    h = 1/n;
    I = speye(n-1,n-1);
    E = sparse(2:n-1,1:n-2,1,n-1,n-1);
    D = E+E'-2*I;
    A = -(kron(D,I)+kron(I,D)) / h^2;
end
