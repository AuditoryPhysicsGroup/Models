function A=schroeder_amp(freq,f1,f2)
% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)

% Calculate amplitude spectrum of a Schroeder harmonic complex
A=zeros(size(freq));
for j=2:numel(freq)
    ff=freq(j);
    if (ff>=f1 && ff<=f2)
        A(j)=1;
    end
end

