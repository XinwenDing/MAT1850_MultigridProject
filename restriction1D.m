function R = restriction1D(n)
    if mod(n,2)~=0
        error('Check n! It should be even.');
    end
    
    R = spdiags(repmat([1, 2, 1], n-1, 1), 0:2, n-2, n);
    R = R(1 : 2 : end, 1:end-1) / 4;
end
