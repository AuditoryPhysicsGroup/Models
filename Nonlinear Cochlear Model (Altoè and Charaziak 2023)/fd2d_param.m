%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428

function param=fd2d_param(dx)

param.dx = dx;               % point spacing in mm (dy = dx)
param.L=5.5;   %the mouse cochlea is 5 mm, we make it a bit longer to suppress low-frequency scattering (see paper)
x=0:dx:param.L;
l=1.85;
fc=80; %natural frequency at the cochlear entrance
param.x=x;
fmin=2;
g=fc*exp(-x/l);
cf=(fc+fmin)*exp(-x/l)-fmin; % we add a "greenwood bend" to the mouse exponential map to suppress low-frequency reflections
param.fr=cf;
param.Nx=length(param.x);
param.d=l;

H=0.33*exp(-x/4/l); % height of mouse cochlea, tapering of height goes with scalae acoustic area
param.H=H;
param.Wscala=param.H; % width of cochlea (same as height, cylindrical model)
[~,i33]=min(abs(x-3));
x33=x(i33);
param.Wbm=[0.095*exp((x(1:i33))/(4*l)) 0.095*exp(x33/(4*l)).*exp((x(i33+1:end)-x33)/(8*l))];%BM width
% param.Wbm=1;
param.rho_fluid=0.001;
[~,i9]=min(abs(cf-10));
c9=cf(i9);
N=sqrt(param.Wbm.*(1./H.^2).*(g)); %Stiffness goes with cf and mass 1/cf
N9=2.75;
param.zweigN=N./N(i9)*N9;  % Parameter N is 2.75 at the ~9 kHz location (based on data), and scales along the cochlea depending on geometry

param.damp=0.25; %damping

Ny=round(param.H/param.dx); %number of vertical points
for i=2:param.Nx-1 %we ensure that at least three consecutive sections have the same number of vertical points, which simplifies finite difference
    if(Ny(i)~=Ny(i-1) && Ny(i)~=Ny(i+1))
       Ny(i)=Ny(i-1); 
    end
end
Ny(end)=Ny(end-1); %we make sure that the last two sections have same height to appropriately apply the boundary conditions in the finite difference
param.Ny=Ny;
