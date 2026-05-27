function param=gerbil_box_model(dx)

param.dx=dx;
[x,cff,mapparams]=cochlear_map(dx,'gerbil');
l=mapparams.d;
param.l=l;
param.x=x;
param.L=x(end);
cf=cff;
param.fr=cf;
param.Nx=length(param.x);
param.d=mapparams.d;


param.Wbm=param.x*0+0.15; % 150 um BM member width

param.rho_fluid=1e-3;% density of water in g/mm3
rho=param.rho_fluid;
zweigN=1.75;

param.H=0.5+param.x*0; %height of the box model chosen to approximate heigth of scala in the base (see Plassmann)
Ny=round(param.H./dx);%not needed, legacy parameter for other things
% [mx,idx]=max(Ny);
% Ny(1:idx)=mx;
for i=2:param.Nx-1
    if(Ny(i)~=Ny(i-1) && Ny(i)~=Ny(i+1))
       Ny(i)=Ny(i-1); 
    end
end
Ny(end)=Ny(end-1);
param.Ny=Ny;
param.Ps=zeros(size(param.x)); %no internal pressure sources
param.Pth=1; 
M=rho*2./param.H;%
mass=param.H*0+M(1)*param.d.^2./(16*zweigN.^2); % calculate mass (constant)

stiffness=mass.*4*pi^2.*cf.^2; 
% stiffness goes with cf^2, this is unrealistic but necessary to get scaling symmetry in a box model 
% this is a gross simplification (see Altoè and Shera 2020 Sci Rep), 
% but perhaps the best choice for the scope of this study 
% as it allows to compare TD solution with frequency-domain solution in a
% simple box (easy to compute, and accurate)


param.M=M;
param.S=param.H.*param.Wbm; % in 2D the area is Wbm*H

param.mass=mass;
param.stiffness=stiffness;
damp_factor=0.2;
param.damp_factor=damp_factor;
param.damping=2*damp_factor*stiffness./(2*pi*cf); %damping coefficient (not factor) of the partition

param.detune_int=1.0; %parameter to detune the second oscillator (vint) in the model
param.int_damp=0.3; % damping coefficient of the second oscillator