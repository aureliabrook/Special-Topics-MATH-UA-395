function x = backward_sub (U, b)
  [m,n] = size(U);
  x = zeros(n,1);        % initialize x to be a column vector
  
  % The last variable calculated directly
  x(n) = b(n)/U(n,n); 
  
  % Back-substitution for remaining variables.
  % ---This could alsobe accomplished without loop by using .*
  for i = (n-1):-1:1 
    x(i) = b(i); % Remember back-sub: x_i = (b_i - sum of A_ij*x_j )/A_ii
    for j = (i+1):n
      x(i) = x(i)-U(i,j)*x(j);
    end
    x(i) = x(i)/U(i,i);
  end
end
