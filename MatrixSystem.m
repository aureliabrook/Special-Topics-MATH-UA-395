function x = MatrixSystem(M,b)
    [A, p] = gauss_eli_srpp(M);
    [n,m] = size(A);
    if n~=m
        error('This function requires a square matrix as an input!')
        return;
    end
    P = zeros(n);
    for i=1:n
        P(i, p(i)) = 1;
    end
    A = P * A;
    U = triu(A);
    L = (tril(A)+diag(ones([n, 1])-diag(tril(A))));
    y = forward_sub(L, P*b);
    x = backward_sub(U, y);
end

