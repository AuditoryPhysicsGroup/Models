function     [Pd0,Vbm,Vcp,Pdbar,alpha,Pd]= fd2d_tapered(params)

  %two-dimensional finite difference solver for tapered geometry 
  %requires at least two consecutive sections with same height
  %works for box models, but in this case it is probably slower than 
  %the solver based on gaussian elimination originally developed by Neely
  %(1981)
  
  Z = params.Z;				% series impedance
  Ybm = params.Ybm;			% BM admittance
  Ycp = params.Ycp;			% Admittance of the entire OoC

  Wbm = params.Wbm;			% width of BM in mm
  H = params.H;		% height of scalae in mm
  dx = params.dx;			% point spacing in mm (dx=dy)
  Nx = params.Nx;			% points in x direction
  Ny = params.Ny;			% points in y direction
  Pth=params.Pth;           % thevei pressure
  Ps=params.Ps.';             % internal pressure sources
  Zch=params.Zch;
  Rstapes=params.Rstapes;
  Rhel=0.01;
%   mult=1i*4*params.omega*dx*params.rho_fluid/Wbm(1);  
  % compute Q vector and Ybm admittances
  % Q = 4*rho*iwo^2*dx * ones([Ny,1]);    % constant stapes displacement
  
  Z0 = Zch(1);
%   Z0=params.Z0;
  % Zstapes=real(Z0*(1+Rstapes)/(1-Rstapes))%+1i*100e-4*params.omega;
  Zstapes=300;
%   real(Z0*(1+Rstapes)/(1-Rstapes))
%   imag(Z0*(1+Rstapes)/(1-Rstapes))/params.omega
  Qin = Pth;%*abs((Zstapes+Z0)/(2*Z0))
  c = Zstapes/(Z(1)*dx);
  Zend=(Zch(end)*(1+Rhel)/(1-Rhel));
%   Yend=Zend;
  cend=Zend/(Z(end)*dx);
  n=Ny(1);
  Q=zeros(sum(Ny),1);
  Q(1:n)= Qin*ones([n,1])*1;    % constant stapes velocity
  p0_idx=zeros(Nx,1);
  %construct the finite difference matrix
  sz=sum(Ny);
  M=sparse(sz,sz);
  kr=1;
  kc=1;
  p0_idx(1)=kc;
  I=speye(n);
  v=-ones(n-2,1);
  u=diag([-2; v],1);
  l=diag([v(1:end); -2],-1);
  A=l+u+4*I;
  A(1,1)=4+(2*Z(1)*H(1)*Ycp(1)*dx);
%    A(1,1)=A(1)+mult*Ybm(1)
  sub=[A*(1+c) -2*I*c];
  M(kc:kc+n-1,kr:kr+2*n-1)=sub*dx;
  kc=kc+n;
  nfuture=Ny(2);
  for i=2:1:(Nx-1)
      npast=n;
      n=nfuture;
      nfuture=Ny(i+1);
      p0_idx(i)=kc;
      Il=-speye(npast);
      Iu=-speye(n);
      v=[-2; -ones(n-2,1)];
      I=speye(n);
      u=diag(v,1);
      l=diag(flipud(v),-1);
      if(n>nfuture)
          Il(end,end)=-2;
          Iu(end,end)=0;
      end
      if(npast>n)
          Il(end,:)=[];
      end
      A=l+u+4*I;
      A(1,1)=A(1)+(2*Z(i).*H(i).*Ycp(i)*dx);
%       A(1,1)=A(1)+mult*Ybm(i);
      sub=[Il A Iu];
      M(kc:kc+size(sub,1)-1,kr:kr+size(sub,2)-1)=sub;
      kr=kr+npast;
      kc=kc+n;
  end
  npast=n;
  p0_idx(end)=kc;
  v=[-2; -ones(n-2,1)];
  I=speye(n);
  u=diag(v,1);
  l=diag(flipud(v),-1);
  Il=2*eye(npast);
  A=l+u+4*I;
  A(1,1)=A(1)+2*(Z(end)*H(end)*Ycp(end)*dx);
%   A(1,1)=A(1)+mult*Ybm(end);
%   sub=[-Il*cend 2*(1+cend)*A];
 sub=[-Il*cend (1+cend)*A];
  M(kc:kc+size(sub,1)-1,kr:kr+size(sub,2)-1)=sub;
% add internal sources
Q(p0_idx)=Q(p0_idx)+Ps;
% solve for pressure
Q=sparse(Q);
% tic;
Pd=M\Q;
% toc;
% calculate quantities
Pd0=Pd(p0_idx);
Vbm = Pd0.*Ybm.'./Wbm.';
Vcp= Pd0.*Ycp.'./Wbm.';
Pdbar=zeros(Nx,1);
for i=1:Nx
    Pdbar(i)=mean(Pd(p0_idx(i):p0_idx(i)+Ny(i)-1));
end
alpha=Pd0./Pdbar;



       
 
 
