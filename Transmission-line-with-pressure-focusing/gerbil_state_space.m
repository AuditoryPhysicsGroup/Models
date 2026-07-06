function [model1D,model2D]=gerbil_state_space(dx,tau)
% dx=spatial resolution in mm
% tau=parameter that controls amplification (amp=0 passive model)
% the function retur the state space representation of the 1D
% transmission-line, and its 2D approximation
param=gerbil_box_model(dx);
wc=2*pi*param.fr;
k0=wc(1)^2*param.M(1);
rme=sqrt(k0*param.mass(1)); %characteristic impedance at entrance (low-freq approximation)
meQ=0.5;
mme=param.M(1);
M=param.M;
model1D.mme=mme;
wca=wc.*param.detune_int; % parameter to "detune" the internal oscillator
nl_knee=1e6*wc(1)./wc; %NL this needs to be calibrated

Ae=sparse(2,2);
Ae(1,1)=-rme./mme;     %first two elements ME state variable
Ae(1,2)=0;
Ae(2,1)=1;
Ap=sparse(4,4);   % for each section, 4 state variable (vbm, ybm, vint and yint)
eff_amp=tau+zeros(param.Nx,1);
m1_loc=round(0.1/dx);
eff_amp(1:m1_loc)=tau*linspace(0,1,m1_loc); % reduce amplification at the cochlear entrance, to prevent standing waves 
% define matrix Ae x'=Ae*x+Be*p
for i=1:param.Nx
    Ap(1,1)=-param.damping(i)./param.mass(i);
    Ap(1,2)=-param.stiffness(i)./param.mass(i);
    Ap(2,1)=1;
    Ap(3,1)=0;
    Ap(3,2)=-eff_amp(i)*wca(i).^2;
    Ap(3,3)=-param.int_damp*wca(i);
    Ap(3,4)=-wca(i)^2;
    Ap(4,3)=1;
    Ae=blkdiag(Ae,Ap);
end
Ap=sparse(1,1);

Ae=blkdiag(Ae,Ap); 

% define matrix Be for a 1D model, with x'=Ae*x+Be*p
Be=sparse(2,1);
Be(1,1)=1./mme;
    Bb=sparse(4,1);
for i=1:param.Nx
    Bb(1,1)=1./param.mass(i);
    Be=blkdiag(Be,Bb);
end
  Bb=sparse(1,1);
  Be=blkdiag(Be,Bb);

% define matrix Ce: acceleration cp=C*x';

Ce=sparse(1,2);
Ce(1,1)=1;
Cc=sparse(1,4);
for i=1:param.Nx
    Cc(1,1)=M(i);
    Cc(1,3)=M(i);
    Ce=blkdiag(Ce,Cc);
end
Cc=sparse(1,1);
Cc(1,1)=1;
Ce=blkdiag(Ce,Cc);

y=zeros(3+4*param.Nx,1); %initial values of pressure and motion
p=zeros(2+param.Nx,1);

Fe=speye(param.Nx+2); %Finite different matrix, N sections plus ME and helicotrema
for i=1%:param.Nx-1
    dSdx=(param.S(i+1)-param.S(i))/param.S(i)/dx;
    Fe(i+1,i)=1./dx^2;
    Fe(i+1,i+1)=-2/dx^2-dSdx/dx/2;
    Fe(i+1,i+2)=1/dx^2+dSdx/dx/2;
end
 %the code is ready to have S depending on distance, functionality not
 %tested in the paper
for i=2:param.Nx-1
    dSdx=(param.S(i+1)-param.S(i))/param.S(i)/dx;
    Fe(i+1,i)=1./dx^2;
    Fe(i+1,i+1)=-2/dx^2-dSdx/dx;
    Fe(i+1,i+2)=1/dx^2+dSdx/dx;
end

for i=param.Nx
    dSdx=0; %no gradients of S at Helicotrema
    Fe(i+1,i)=1./dx^2-dSdx/dx;
    Fe(i+1,i+1)=-2/dx^2;
    Fe(i+1,i+2)=1/dx^2+dSdx/dx;
end

Fe(param.Nx+2,param.Nx+2)=1;
Fe(1,1)=-1./dx;
Fe(1,2)=1./dx;
model1D.q=sparse(param.Nx+2,1);
Fmin1=inv(Fe);
Common=inv(speye(numel(y))-Be*Fmin1*Ce);

model1D.A=(Common*Ae);
model1D.B=(Common*Be);
model1D.Ae=Ae;
model1D.Be=Be;
model1D.Ce=Ce;
model1D.Fcb=Fe-Ce*Be;
model1D.FcbFast=decomposition(model1D.Fcb,'lu'); 
%in our system matlab solver is now faster solving Fy=x, with tridiagonal
%solver than LU decomposition
model1D.dx=dx;
model1D.Nx=param.Nx;
model1D.ystart=y;
model1D.pstart=p;
model1D.cf=param.fr;


%% 2D approximation, all matrixes the same,except BE
model2D=model1D;
%
% Parameter to ensure long-wave behavior near entrance and helicotrema
% (first 8 sections are 1D)
skip_h=8; 
H=param.H;

% define matrixes D and E such that E*p0=D*pbar;
% A1 and B1 
D1=speye(param.Nx+2,param.Nx+2);

Nx=param.Nx;
Nx=param.Nx;
hx2=H.^2./dx^2*4/9;
hx2(1:skip_h)=0; hx2(end-skip_h:end)=0;
hx2=[0 hx2 0].';

D22=spdiags([-circshift(hx2,-1) 2*hx2 -circshift(hx2,1)],-1:1,Nx+2,Nx+2); %2nd spatial derivative
hx4=H.^4./dx^4*1/63;
hx4(1:skip_h)=0; hx4(end-skip_h:end)=0;
hx4=[0 hx4 0].';
D44=spdiags([circshift(hx4,-2) -4*circshift(hx4,-1) ...
    6*hx4 -4*circshift(hx4,1) circshift(hx4,2)],-2:2,Nx+2,Nx+2); %4ht spatial derivative

D=D1+D22+D44; % sum them all

E1=speye(Nx+2,Nx+2);
hx2=H.^2./dx^2*1/9; 
hx2(1:skip_h)=0; hx2(end-skip_h:end)=0;
hx2=[0 hx2 0].';
E22=spdiags([-circshift(hx2,-1) 2*hx2 -circshift(hx2,1)],-1:1,Nx+2,Nx+2);
hx4=H.^4./dx^4*1/945;
hx4(1:skip_h)=0; hx4(end-skip_h:end)=0;

hx4=[0 hx4 0].';
E44=spdiags([circshift(hx4,-2) -4*circshift(hx4,-1) ...
    6*hx4 -4*circshift(hx4,1) circshift(hx4,2)],-2:2,Nx+2,Nx+2);
E=E1+E22+E44;
%Be2d=Be1*E^(-1)*D with S matrix representation of the "pressure focusing" spatial filter 
% that relates driving and average pressure, S=A1^(-1)B1. See paper
Be2d=Be*(E\D);

Common=inv(speye(numel(y))-Be2d*Fmin1*Ce);
model2D.A=Common*Ae;
model2D.B=Common*Be2d;
model2D.dx=dx;
model2D.Nx=param.Nx;

N=numel(Fe(1,:)); 
model2D.N=N;

model2D.Ae=Ae;
model2D.Be=Be2d;
model2D.Ce=Ce;

model2D.E=E;
model2D.D=D;
model2D.Be1d=Be;

model2D.Ee=decomposition(model2D.E,'lu'); % save some time using lu decomposition
LargeMatrix=([Fe -Ce*Be; -D E]);
% this is to solve (Fe-CBE^(-1)D)^(-1)y=x avoiding matrix inversion 
% the solution is to solve Fy-CBu=x -Ey+Du=0
% numerically we solve LargeMatrix * y2=[x 0]^T
% and then discard half of the solution
p=reshape([1:N; N+1:2*N], 1, []);
LargeReshaped=LargeMatrix(p,p);
% we reshape the matrix so it is banded
model2D.p=p;
model2D.bigD=sparse(N,2*N);
model2D.bigD(:,1:2:end)=model2D.D;
% this big D is like D but takes only the part of the solution we need.
[model2D.L,model2D.U,model2D.P] = lu(LargeReshaped);
% use LU decomposition to speed up matters
% would be great to have a low-level solver for this banded matrix
[~, pLU] = max(model2D.P, [], 2);
pp2=p(pLU);
model2D.pp2=pp2;






