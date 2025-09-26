function param=param_solver(param_in,f,tau)

param=param_in;
param.tau=tau+param.fr*0;
ramp=10;
taper=linspace(1,-0,ramp);
taperup=linspace(0,1,ramp);
param.tau(1:ramp)=param.tau(1:ramp).*taperup;

param.tau(end-ramp+1:end)=param.tau(end-ramp+1:end).*taper;
param.omega=2*pi*f;
param.Rstapes=0.1;
param.Z=Zseries(param,f);
[ZYcp,ZYbm]=ZseriesYcp(param,f);
param.Ycp=ZYcp./param.Z;
param.Ybm=ZYbm./param.Z;
param.Zch = csqrt([param.Z(1)/param.Ycp(1),param.Z(end)/param.Ycp(end)]);
% param.Z0=sqrt(param.Z(1)/param.Ybm(1));