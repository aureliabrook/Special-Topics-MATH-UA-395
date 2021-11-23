function y = forward_sub (L, b)
  [m,n] = size(L);
  y = zeros(n,1);        % initialize x to be a column vector
  
  % The first variable calculated directly
  y(1) = b(1)/L(1,1); 
  
  % Back-substitution for remaining variables.
  % ---This could also be accomplished without loop by using .*
  for i = 1:n
    y(i) = b(i); % Remember forward-sub: y_i = (b_i - sum of A_ij*x_j )/A_ii
    for j = 1:i-1
      y(i) = y(i)-L(i,j)*y(j);
    end
    y(i) = y(i)/L(i,i);
  end

end