function Z = fitKnownSines(xi,yi,f0,persist)
% function Z = fitKnownSines(xi,yi,f0,?persist=false?)
% Fits data xi,yi to sinusoids of known frequency f0, where
%   f0 can be an array.
% Returns an array of complex numbers.
% Sinusoid amplitudes are abs(Z), phases are angle(Z)
% If persist then saves/reuses the matrix M
  
  persistent M
  
  if nargin<4
    persist = false;
  end
  
  if isempty(M) || ~persist
    phi = 2*pi*repmat(xi(:)-xi(1),1,2*numel(f0)) .* repelem(torow(f0),2);
    M = zeros(size(phi));
  
    M(:,1:2:end) =  cos(phi(:,1:2:end));
    M(:,2:2:end) = -sin(phi(:,2:2:end));
  end

  ri = M\yi(:);
      
  Z = reshape(complex(ri(1:2:end),ri(2:2:end)),size(f0));
  
  if ~persist
    M = [];
  end
    
end