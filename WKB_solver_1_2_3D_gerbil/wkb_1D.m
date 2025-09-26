%This code accompanies the manuscript by  Altoè A. "The role of scalae geometry in cochlear mechanical responses: pressure focusing in two and three dimensions."
%under review
function     [k,Vbm,Vcp,Pdbar,alpha,Pd]= wkb_1D(params)

  %1D WKB solver

  Z = params.Z;				% series impedance
  Ybm = params.Ybm;			% BM admittance
  Ycp = params.Ycp;			% Admittance of the entire OoC
  
  Wbm = params.Wbm;			% width of BM in mm
  dx = params.dx;			% point spacing in mm (dx=dy)
  Pth=params.Pth;           % thevei pressure

  S=params.S;

  wkb_fac=248.3/2; %empirical factor to get same amplitudes as with finite differences
  
  ZY=Z.*Ycp;
  k2lw=(-ZY);


  alpha=1;
 
  k = csqrt(alpha.*k2lw);
  int=cumtrapz(-1i*k*dx);
  Pdbar=wkb_fac.*Pth.*sqrt((S(1).*k(1))./(S.*k)).*exp(int); %average pressure
  Pd=Pdbar.*alpha; %pressure across partition
  Vcp=Pd.*Ycp; %Partition's velocity
  Vbm=Pd.*Ybm./Wbm; %BM velocity 
  Pd0=Pd;
 
 
