function y = csqrt(z)
% function y = csqrt(z) -- sqrt that avoids branch cuts when
%   called with complex array z. Uses clog(z).
%  
% C.A.Shera
  
  y = exp(clog(z)/2);
  
  return
  
function y = clog(z)
% function y = clog(z) -- log that uses the unwrapped phase to
%   avoids branch cuts in the complex plane 
%
% C.A.Shera

  logz = log(z);
  y = real(logz) + 1i*unwrap(imag(logz));
  
  return
  