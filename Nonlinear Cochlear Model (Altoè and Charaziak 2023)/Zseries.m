%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428

function z = Zseries(param,f)
% function z = Zseries
% The values of unspecified arguments are
% taken from corresponding global variables

S=param.H.^2*pi/2;

 z = 1i*2*pi*f*(param.rho_fluid)./(S);


