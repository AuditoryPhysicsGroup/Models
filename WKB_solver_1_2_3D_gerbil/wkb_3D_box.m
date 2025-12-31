%This code accompanies the manuscript by  Altoè A. "The role of scalae geometry in cochlear mechanical responses: pressure focusing in two and three dimensions."
%under review
function     [k,Vbm,Vcp,Pdbar,alpha,Pd]= wkb_3D_box(params)

  %3D WKB solver for a box model (with acoustic dimensions matched to a
  %cylindrical model)
  Z = params.Z;				% series impedance
  Ybm = params.Ybm;			% BM admittance
  Ycp = params.Ycp;			% Admittance of the entire OoC
  
  area=((params.R.^2-(params.Wbm/pi).^2)*pi); %area of the duct of the cylindrical model 
  W=sqrt(area); %we keep the same area for the box
  H=W/2;
  dx = params.dx;			% point spacing in mm (dx=dy)

  Pth=params.Pth;           % thevei pressure
 
  S=params.S;
  wkb_fac=1; %amplitude scaling factor (not used)
  
  ZY=Z.*Ycp;
  k2lw=(-ZY);
  klw=csqrt(k2lw);
  b=W; %use de Boer notation
  Wbm=params.Wbm*(pi/2);
  e=(Wbm)./W;
  dw=0.5; %partition centered along the y-axis



  alpha=klw*0+1; %start from alpha=1 everywhere
    for n=1:10 %usual iteration for computing alpha [see Altoè & Shera (2020) Sci. Rep.]
    k= csqrt(alpha.*k2lw);
    k = abs(real(k)) + 1i*abs(imag(k)).*sign(imag(klw));
    m0=k;
    % m02 = abs(real(m02)) + 1i*abs(imag(k)).*sign(imag(klw));
    c0=2*e/pi;
    alpha=k.*H./tanh(k.*H);
    m02=alpha.*k2lw;
    av0=(m02.*tanh(m0.*H)./c0);
    for i=1:30 % 30 modes should do
        l=pi*i./b;
        m=csqrt(k.^2+l.^2);
        m = abs(real(m))+1i*abs(imag(m)).*sign(imag(klw));
        cn=-4*e.*cos(i*e*pi/2).*cos(pi*i*dw)./(e.^2.*i^2-1)/pi;
        alpha=alpha+av0.*cn./(m.*tanh(m.*H)).*sin(pi*i*e/2).*cos(pi*i*dw)./(pi.*i*e/2);
    end
    end
  
  int=cumsum(-1i*k*dx);
  Pdbar=wkb_fac.*Pth.*sqrt((S(1).*k(1))./(S.*k)).*exp(int); %average pressure
  Pd=Pdbar.*alpha; %pressure across partition
  Vcp=Pd.*Ycp; %Partition's velocity
  Vbm=Pd.*Ybm; %BM velocity 

 
 
