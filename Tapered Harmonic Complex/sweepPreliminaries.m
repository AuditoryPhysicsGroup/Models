function [freqs,R,T]=sweepPreliminaries(r,fs)
% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)

% preliminary calculation to design a periodic sweep/Schroeder's complex

Nf=round(fs/r);
dt=1./fs;
T=Nf*dt;
R=1./T;
freqs=0:R:fs/2;
