%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428

function g=me_gain(f)

freq=f;
fme=15;
b=freq./fme;
g=1i*b./(1+1i*b).*exp(-1i*2*pi*freq*40e-3)/2;
