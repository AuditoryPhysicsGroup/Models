function [zycp,zybm] = ZseriesYcp2D(param,f)
% function zy = zseriesybm
% The values of unspecified arguments are
% taken from corresponding global variables

b = f./param.fr;
sc=10;

% bb = b.*b;
% fourNd2=(param.zweigN.*4./param.d).^2;
% t0=param.tau;
% zybm = - (fourNd2.*bb).*(1)./(1-bb+1i*2*param.damp.*b); %mouse 1.05 live
% zycp= zybm.*(1+1i*b.*t0);
om=2*pi*f;
om2=om^2;
% zybm=-(param.M.*param.Wbm).*om^2./(-param.mass.*om2+1i*param.damping.*om+param.stiffness);
zybm=-(param.M).*om^2./(-param.mass.*om2+1i*param.damping.*om+param.stiffness);

% param.detune_int=1.;
% param.int_damp=0.5;
om_amp=2*pi*param.fr*param.detune_int;
% zycp=zybm.*(1+sc*1i*(param.tau.*bfve)./(-bfve.^2+1.5*1i*bfve+1));
damp=param.int_damp*om_amp;
% zycp=zybm.*(1+(param.detune_int*1i*om*param.tau)./(om_amp+1i*om));%./(-om.^2+1i*damp*om.*om_amp+om_amp.^2));
% zycp=zybm.*(1+(param.detune_int*om_amp.*1i.*om.*param.tau)./(om_amp.^2+1i*damp-om.^2));
% this works well with detune_int large (10) and damp large (1)
zycp=zybm.*(1-(om_amp.^2.*param.tau*(1-0.0*1i))./(om_amp.^2+1i*damp.*om-om.^2));%./(-om.^2+1i*damp*om.*om_amp+om_amp.^2));

% zybm = - (fourNd2.*bb).*(1-t0*0.15)./(1-bb+1i*2*param.damp.*b); %mouse 1.05 live
