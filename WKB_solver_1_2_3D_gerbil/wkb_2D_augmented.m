%This code accompanies the manuscript by  Altoè and Charaziak (2023) "Characterizing Nonlinear Wave Amplification in the Mouse Cochlea" 

function     [k,Vbm,Vcp,Pdbar,alpha,Pd]= wkb_2D_augmented(params)

  %two-dimensional WKB solver for a 2D model with same fluid area of a 3D
  %model [see paper]
  
  Z = params.Z;				% series impedance
  Ybm = params.Ybm;			% BM admittance
  Ycp = params.Ycp;			% Admittance of the entire OoC
  
  dx = params.dx;			% point spacing in mm (dx=dy)

  Pth=params.Pth;           % input pressure
  
  S=params.S;
 
  wkb_fac=1; %empirical scaling factor (not used)
  
  ZY=Z.*Ycp;
  k2lw=(-ZY);
  klw=csqrt(k2lw);
  area=((params.R.^2-(params.Wbm/pi).^2)*pi/2); %area of the duct % height of scalae in mm %3D correction
  Wbm=params.Wbm*(pi/2);
  H=area./(Wbm);
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
  Vbm=Pd.*Ybm; %BM velocity 

 
 
