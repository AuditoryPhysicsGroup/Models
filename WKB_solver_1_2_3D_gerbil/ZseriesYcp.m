function [zycp,zybm] = ZseriesYcp(param,f)
% function zy = zseriesybm
% The values of unspecified arguments are
% taken from corresponding global variables

b = f./param.fr;
% bb = b.*b;
% fourNd2=(param.zweigN.*4./param.d).^2;
% t0=param.tau;
% zybm = - (fourNd2.*bb).*(1)./(1-bb+1i*2*param.damp.*b); %mouse 1.05 live
% zycp= zybm.*(1+1i*b.*t0);
om=2*pi*f;
om2=om^2;
zybm=-(param.M.*param.Wbm).*om^2./(-param.mass.*om2+1i*param.damping.*om+param.stiffness);
zycp=zybm.*(1+param.tau);
% zybm = - (fourNd2.*bb).*(1-t0*0.15)./(1-bb+1i*2*param.damp.*b); %mouse 1.05 live
