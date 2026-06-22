function [spec,phase]=schroeder_phase_complex_spectrum(A,C)
% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)



%Calculate phase from spectrum according to Schroeder's equations
%M. Schroeder IEEE trans. Inf. Theory 16, 85–89 (1970)
%but with the paramter C to control the crest factor and phase curvature
%for a spectrum with uniform amplitude within f1 and f2, this script
%generates the "Schroeder-phase harmonic" detailed by
% Lentz, Jennifer J., and Marjorie R. Leek (2001) 
% Journal of the Association for Research in Otolaryngology 2.4 (2001): 408-422.
power=abs(A.^2);
pk=power./sum(power);
tk=cumsum(pk);
phase=2*pi*C*cumsum(tk);
spec=abs(A).*exp(1i*(phase(1:end)));

