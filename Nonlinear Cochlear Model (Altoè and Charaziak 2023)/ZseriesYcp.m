%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428

function [zycp,zybm] = ZseriesYcp(param,f)
%calculate the product of Z*Y for BM and partition admittance

b = f./param.fr;
bb = b.*b;
fourNd2=(param.zweigN.*4./param.d).^2;
gdy=param.gdy;
zybm = - (fourNd2.*bb).*(1)./(1-bb+1i*2*param.damp.*b); %mouse 1.05 live
zycp= zybm.*(1+gdy);
