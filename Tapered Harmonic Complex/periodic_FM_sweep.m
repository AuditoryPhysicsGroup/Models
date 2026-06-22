function [sig,dfdtR]=periodic_FM_sweep(f1,f2,rate,fs)
% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)

% design one period of a periodic FM sweep, making sure that the phase is
% continuous is at the beginning/end of each cycle
  dt=1/fs;
  T=1/rate;
  dfdt=(f2-f1)/T;
  T=(T);
  phiCend=f1*T+dfdt*T^2/2;
  c=ceil(phiCend);
  dfdtR=2*(c-f1*T)/T^2;
  t=0:dt:(T-dt);
  t=t-floor(t/T)*T;
  phi=2*pi*(f1*t+dfdtR*t.^2/2);
  phi(end)
  sig=sin(phi);

  
  