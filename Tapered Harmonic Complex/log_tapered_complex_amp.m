function A=log_tapered_complex_amp(freq,f1,f2,order)
% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)


% Calculate the amplitude of a Schroeder complex with frequency that
% increases exponentially over time
% I.e., periodic version of log-sweep
A=zeros(size(freq));
A(2:end)=sqrt(f1./freq(2:end));
for j=2:numel(freq)-1
     ff=freq(j);
    if ff<f1
        A(j)=A(j)*(ff/f1)^order; 
    elseif ff>f2
        A(j)=A(j)*(f2/ff)^order;
   % else
    %    A(j)=1;
    end
end
