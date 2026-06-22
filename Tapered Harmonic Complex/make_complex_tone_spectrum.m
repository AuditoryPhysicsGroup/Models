function [spec,phase]=make_complex_tone_spectrum(A)
% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)



%Calculate phase from spectrum according to Schroeder's equations
%M. Schroeder IEEE trans. Inf. Theory 16, 85–89 (1970)
power=abs(A.^2);
pk=power./sum(power);
tk=cumsum(pk);
phase=-2*pi*cumsum(tk);
spec=abs(A).*exp(1i*(phase(1:end)));

