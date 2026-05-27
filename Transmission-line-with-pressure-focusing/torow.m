function y = torow(x)
% function y = torow(x)
% return x as a row vector
  
  [m,n] = size(x);
  if (~isempty(x) && (m ~= 1 & n ~= 1))
    error ('x not a vector');
  end
  
  if (n>=m)
    y = x;
  else
    y = x.';
  end
  
  return
  