% Matlab script to calculate cochlear model responses
% presented in the paper 
%    Altoe A, Shera CA (2020) The cochlear ear horn: Geometric origin of
%    tonotopic variations in auditory signal processing. Sci Rep 10:20528
%    https://doi.org/10.1038/s41598-020-77042-w
% 

clear all;
close all;
clc;

L=22.5;                       % length of cat's cochlea
dx=L/1150;                    % 1150 sections
x=dx:dx:L;
dx=x(2)-(x(1));
l=L/2.1/log(10);              % greenwood L to Zweig l
fb=57.4;                      % higher CF at the entrance
pos_sc=650;                   % apical/basal transition location
cff=fb.*10^(-2.1)*(10.^((L-x)./L*2.1)-0.8); % greenwood map
cf=cff(1)*exp(-x/l);          % exponential frequency-map used to compute parameters 
h=L*((cf./cf(1)).^(2/3)*0.047+0.015); % height variaton
H=h(1);                       % height at the entrance
rho=1.;                       % fluid density
Mb=rho/H.^2/2;                % acoustic fluid mass at the entrance
Nst=1*pi/2;                   % N=pi/2*tau_wf in the base
N=x*0+Nst;                    % make N a vector 
N(pos_sc:end)=Nst.*exp(-(x(pos_sc:end)-x(round(pos_sc)))/4/l); % delay scaling break in the apex
sp_gain=exp(x/2/l);           % pressure gain due to area tapering in the base
sp_gain(pos_sc:end)=exp(((x(pos_sc:end)-x(pos_sc))/4/l)).*sp_gain((pos_sc)); % and in the apex
sp_gain=sp_gain';
mratio1=16*N.^2./l.^2;        % fluid/CP mass ratio at the cochlear entrance computed from the delay
mb=Mb(1)./mratio1(1);         % CP mass at the entrance
st=(cf(1)).^2.*mb.*cf./cf(1); % CP stiffness along the cochlea (exponential)

Ma=Mb./(cf/cf(1));            % acoustic mass of fluids tapering (e^x/l in base)
Ma(pos_sc:end)=Ma(pos_sc).*sqrt(cf(pos_sc)./cf(pos_sc:end)); % (e^x/2l in the apex)
ma=st./(cff).^2;              % correct acoustic mass in the apex for greenwood map's "bend"
mratio=(Ma./ma);              % mass ratio along the cochlea
cf=cff;                       % replace the exponential frequency map with the true one

tau=1.3;                      % parameter that controls amplification (tau=0 for passive)
d=2*0.12;                     % damping
h=h';

itt=0;
fflist=logspace(log10(0.01),log10(50),1e3);

tic;
for ff=fflist
  itt=itt+1;
  beta=ff./cf;                % normalized frequency
  % numerator of long-wave wavenumber^2...
  k_2_n=mratio.*beta.^2.*(1+tau.*beta.*1i); 
  % denominator of long-wave wavenumber^2...
  k_2_d=(1-beta.^2+1i.*d.*beta); 
  % CP admittance
  admittance=1i*beta.*(1+tau.*beta.*1i)./(1-beta.^2+1i.*d.*beta); 
  k_2=(k_2_n./k_2_d).';
  k_lw=csqrt(k_2);            % long wave approximation
  k=k_lw;                     % start from long-wave approximation 
  alpha_WKB=k.*h./(tanh(k.*h));

  % 20-30 iterations are sufficient to get a good approximation
  for i=1:25                            
    alpha_WKB=((k).*h)./tanh(k.*h);
    kp=k;
    k=csqrt((k_2).*alpha_WKB);
    % force the algorithm to pick the solution 
    % of a forward traveling wave with imaginary wavenumber sign as that of
    % the long-wave approximation (strongly reduces artifacts above CF)
    k=abs(real(k))+1i.*abs(imag(k)).*sign(imag(k_lw)); 
  end

  integral=cumtrapz(k.*dx);
  Pbar=sp_gain.*csqrt(k(1)./k).*exp(-1i*integral); % averaged pressure
  P0=(alpha_WKB).*Pbar;       % driving pressure
  BM=P0.*admittance.';        % BM response

  alpha_sol(:,itt)=alpha_WKB(1:1:end); % save the solutions (location,frequency)
  BM_mag(:,itt)=abs(BM(1:1:end));
  BM_angle(:,itt)=unwrap(angle(BM));
  Pd(:,itt)=P0;
  ksol(:,itt)=k;
end
toc;

probe_locations=[200:100:1.1e3];% plot response for these locations
figure(1);
for i=probe_locations
  BM_gain=mag2db((BM_mag(i,:)));
  maxGain=max(BM_gain);
  % plot only for frequency where the response is within 80 dB of peak 
  % (for clarity)
  okk=find(BM_gain>(maxGain-80)); 
  subplot(2,2,1);
  semilogx(fflist(okk),BM_gain(okk)); hold all 
  subplot(2,2,2);
  semilogx(fflist(okk),(BM_angle(i,okk))/2/pi); hold all 
  subplot(2,2,3);
  semilogx(fflist(okk),mag2db((abs(alpha_sol(i,okk))))); hold all 
  subplot(2,2,4);
  semilogx(fflist(okk),unwrap(angle(alpha_sol(i,okk)))/2/pi); hold all 

end

subplot(2,2,1);
xlabel('Frequency [kHz]');
ylabel('BM gain [dB re e.c. pressure]');
xlim([.1 30]);
ylim([-10 50]);
set(gca,'FontSize',12);

subplot(2,2,2);
xlabel('Frequency [kHz]');
ylabel('BM phase [cycles]');
xlim([.1 30]);
ylim([-3 0.5]);
set(gca,'FontSize',12);

subplot(2,2,3);
xlabel('Frequency [kHz]');
ylabel('|alpha| [dB]');
xlim([.1 30]);
ylim([-5 30]);
set(gca,'FontSize',12);

subplot(2,2,4);
xlabel('Frequency [kHz]');
ylabel('\angle{alpha} [cycles]');
xlim([.1 25]);
ylim([-0.25 0.1]);
set(gca,'FontSize',12);



