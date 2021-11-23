% Perform gaussian elimination.
% In place of changing elements of A to zero. Use them to save entries of L
% Returns the permutation vector and matrix A after transformation.

function [A, p] = gauss_eli_srpp (A)
  [n,m] = size(A);
  if n~=m
    error('This function requires a square matrix as an input!')
    return;
  end
  p = (1:n)';	 % the permutation vector
  s = max(abs(A'));    % scale for each row
  
  for k = 1:(n-1) % k-th step of Gaussian Elimination
    d_s = abs(A(p(k),k))/s(p(k));
    pivot_index = k;
    
    % Find next Pivot row
    for i = k+1:n
      z = abs(A(p(i),k))/s(p(i)); 
      if z > d_s,  
        d_s = z;  
        pivot_index = i;  
      end
    end    
    
    % Change the permutation vector indicating row-swaps
    curr_index = p(k); 
    p(k)= p(pivot_index);  
    p(pivot_index) = curr_index; 
    
    % Gaussian Elimination for step-k. 
    for i = (k+1):n
      m = A(p(i),k)/A(p(k),k); % The multiplier for row k+1, ...
      
      %change elements of every column in row i
      for j = (k+1):n
        A(p(i),j) = A(p(i),j)- m*A(p(k),j);
      end
      
      A(p(i),k) = m; % We can save the multiplier right in A in place of L
    end
    
  end
end
