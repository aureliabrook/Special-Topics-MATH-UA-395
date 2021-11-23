function I = MatrixInverse(M)
    [n, m] = size(M);
    I = zeros(n);
    for i = 1:n
        curr = zeros([n, 1]);
        curr(i) = 1;
        b=MatrixSystem(M, curr);
        for j=1:n
            I(j, i) = b(j);
        end
    end
    disp(M*I);
end

