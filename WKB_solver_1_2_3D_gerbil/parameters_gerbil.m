function param=parameters_gerbil(dx)

param.dx=dx;
[x,cff,mapparams]=cochlear_map(dx,'gerbil');
l=mapparams.d;
param.x=x;
param.L=x(end);
cf=cff;
param.fr=cf;
param.Nx=length(param.x);
param.d=mapparams.d;

load BM_pol.mat
b=polyval(pp,param.x)*1e-3;%*2/pi;
[~,ab]=min(abs(x-0.4*param.L));
xab=param.x(ab);

param.ab=ab;
param.Wbm=b; 
param.rho_fluid=1e-3;% density of water in g/mm3
rho=param.rho_fluid;
r0=b/pi;
S=polyval(sp,param.x);
param.R=sqrt(4*S/pi);
Seff=pi*(param.R.^2-r0.^2)/4; %effective area


param.Ny=param.x*0+1;%not needed, legacy parameter for other things

param.Ps=zeros(size(param.x)); %no internal sources
param.Pth=1; 
init_st=1.5e9*1e-3;%1.5e9 Pa/m in the 40 kHz region (Dong Olson 2009)
stiffness=init_st.*cf./cf(1); %scale it along the cochlea
mass=stiffness./(2*pi*cf).^2; % CP mass from resonator theory


M=rho*param.Wbm./Seff;%
param.M=M;
param.S=S;

param.mass=mass;
param.stiffness=stiffness;
damp_factor=0.21;
param.damp_factor=damp_factor;
param.damping=2*damp_factor*stiffness./(2*pi*cf); %damping coefficient (not factor)