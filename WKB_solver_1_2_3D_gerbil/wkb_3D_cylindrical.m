%This code accompanies the manuscript by  Altoè A. "The role of scalae geometry in cochlear mechanical responses: pressure focusing in two and three dimensions."
%under review
function     [k,Vbm,Vcp,Pdbar,alpha,Pd]= wkb_3D_cylindrical(params)

  %3D WKB solver for a cylindrica model
  
  Z = params.Z;				% series impedance
  Ybm = params.Ybm;			% BM admittance
  Ycp = params.Ycp;			% Admittance of the entire OoC
  Wbm = params.Wbm;			% width of BM in mm
  dx = params.dx;			% point spacing in mm (dx=dy)

  Pth=params.Pth;           % thevei pressure

  S=params.S;
  wkb_fac=1; %empirical scaling factor % not used
  ZY=Z.*Ycp;
  k2lw=(-ZY);
  klw=csqrt(k2lw);
  r1=params.Wbm/pi;
  r2=params.R;

  area=pi*r2.^2-pi*r1.^2;
  alpha=klw*0+1; %start from alpha=1 everywhere
    for n=1:30 %usual iteration for computing alpha [see Altoè & Shera (2020) Sci. Rep.]
    k= csqrt(alpha.*k2lw);
        k = abs(real(k)) + 1i*abs(imag(k)).*sign(imag(klw));
    pD0=besselk(0,k.*r1)+besselk(1,r2.*k)./besseli(1,r2.*k).*besseli(0,k.*r1);
    pDbar=pi*2*(r1.*besselk(1,k.*r1)-r2.*besselk(1,k.*r2)+...
        besselk(1,r2.*k)./besseli(1,r2.*k).*(-r1.*besseli(1,k.*r1)+r2.*besseli(1,k.*r2)))./k./area;
    alpha=pD0./pDbar;
    end
   
  int=cumsum(-1i*k*dx);
  Pdbar=wkb_fac.*Pth.*sqrt((S(1).*k(1))./(S.*k)).*exp(int); %average pressure
  Pd=Pdbar.*alpha; %pressure across partition
  Vcp=Pd.*Ycp; %Partition's velocity
  Vbm=Pd.*Ybm./Wbm; %BM velocity 

 
 
