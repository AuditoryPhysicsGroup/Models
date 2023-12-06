%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428

function param=param_solver(param_in,f,gdy)
% prepares parameters for calculatuing the model's solution
% the param_in are the geometrical paramters (from fd2d_param)
% gdy can be a vector or a scalar, and represent  g*dy as in Eq.(20) in the
% paper
param=param_in;
param.gdy=gdy+zeros(size(param.x));
param.Pth=0;
param.Ps=zeros(size(param.x)); %internal pressure sources, initialized at zero

end_sp=round(0.1/param.dx); % kill activity over the last 0.1 mm tofurther prevent scattering from helicotrema using fd2d
ramp=round(end_sp-1);
taper=linspace(1,0,ramp);
taper=[taper zeros(1,end_sp-ramp)];
param.gdy(end-end_sp+1:end)=param.gdy(end-end_sp+1:end).*taper;
param.omega=2*pi*f;
param.Rstapes=0.;
param.Z=Zseries(param,f);
[ZYcp,ZYbm]=ZseriesYcp(param,f);
param.Ycp=ZYcp./param.Z;
param.Ybm=ZYbm./param.Z;
param.Zch = sqrt([param.Z(1)/param.Ycp(1),param.Z(end)/param.Ycp(end)]);