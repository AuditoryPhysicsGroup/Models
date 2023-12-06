%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428

function     [Pd0,Vbm,Vcp,Pdbar,alpha,Pd]= wkb_tapered(params)

  %two-dimensional WKB solver based on the work of Duifhuis (1988)
  
  Z = params.Z;				% series impedance
  Ybm = params.Ybm;			% BM admittance
  Ycp = params.Ycp;			% Admittance of the entire OoC
  
  Wbm = params.Wbm;			% width of BM in mm
  H = params.H;				% height of scalae in mm
  dx = params.dx;			% point spacing in mm (dx=dy)
  Nx = params.Nx;			% points in x direction
  Ny = params.Ny;			% points in y direction
  Pth=params.Pth;           % thevei pressure
  Ps=params.Ps.';             % internal pressure sources
  Zch=params.Zch;
  Rstapes=params.Rstapes;
  Rhel=0.01;
  S=pi.*H.^2;
  Pd0=0;
  wkb_fac=248.3/2; %empirical factor to get same amplitudes as with finite differences

  ZY=Z.*Ycp;
  k2lw=(-ZY);
  klw=sqrt(k2lw);

  alpha=klw*0+1; %start from alpha=1 everywhere
    for n=1:8 %usual iteration for computing alpha [see Altoè & Shera (2020) Sci. Rep.]
    kH = sqrt(alpha.*k2lw).*H;
    kH = abs(real(kH)) + 1i*abs(imag(kH)).*sign(imag(klw));
    alpha = real(kH)./tanh(real(kH));
    end
    % two additional iterations based on series expansion of alpha (see Shera & Altoe 2023, PNAS)
    % to reduce numerical noise in cut-off region
    for n=1:2 
    kH = sqrt(alpha.*k2lw).*H;
    kHr = real(kH); kHi = imag(kH);
    % series expansion in kHi about kHi=0... 
    ar = kHr.*coth(kHr)-(kHr.*coth(kHr)-1).*csch(kHr).^2.*kHi.^2;
    ai = (coth(kHr)-kHr.*csch(kHr).^2).*kHi + ...
         (kHi.^3/6).*csch(kHr).^4.*(4*kHr+2*kHr.*cosh(2*kHr)-3*sinh(2*kHr));
    alpha = ar + 1i*ai;
    k = sqrt(alpha.*k2lw);
    end
  int=cumtrapz(-1i*k*dx);
  Pdbar=wkb_fac.*Pth.*sqrt((S(1).*k(1))./(S.*k)).*exp(int); %average pressure
  Pd=Pdbar.*alpha; %pressure across partition
  Vcp=Pd.*Ycp; %Partition's velocity
  Vbm=Pd.*Ybm; %BM velocity 
       
 
 
