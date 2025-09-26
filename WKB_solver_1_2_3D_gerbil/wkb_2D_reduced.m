%This code accompanies the manuscript by  Altoè A. "The role of scalae geometry in cochlear mechanical responses: pressure focusing in two and three dimensions."
%under review
function     [k,Vbm,Vcp,Pdbar,alpha,Pd]= wkb_2D_reduced(params)

  %two-dimensional WKB solver, for a 2D model with same "height" of a 3D
  %model
  Z = params.Z;				% series impedance
  Ybm = params.Ybm;			% BM admittance
  Ycp = params.Ycp;			% Admittance of the entire OoC
  H=params.R;
  Wbm = params.Wbm;			% width of BM in mm
  dx = params.dx;			% point spacing in mm (dx=dy)

  Pth=params.Pth;           % thevei pressure
  
  S=params.S;

  wkb_fac=1; %empirical scaling factor (not used)
  ZY=Z.*Ycp;
  k2lw=(-ZY);
  klw=csqrt(k2lw);

  alpha=klw*0+1; %start from alpha=1 everywhere
    for n=1:30 %usual iteration for computing alpha [see Altoè & Shera (2020) Sci. Rep.]
    k= csqrt(alpha.*k2lw);
    k = abs(real(k)) + 1i*abs(imag(k)).*sign(imag(klw));
    alpha = (k.*H)./tanh((k.*H));
    end
  int=cumtrapz(-1i*k*dx);
  Pdbar=wkb_fac.*Pth.*sqrt((S(1).*k(1))./(S.*k)).*exp(int); %average pressure
  Pd=Pdbar.*alpha; %pressure across partition
  Vcp=Pd.*Ycp; %Partition's velocity
  Vbm=Pd.*Ybm./Wbm; %BM velocity 

 
