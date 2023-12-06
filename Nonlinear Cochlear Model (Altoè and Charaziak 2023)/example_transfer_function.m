%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428

%script to calculate BM and TM transfer functions in the model.
clear all;
close all;
clc;
fac=1;
dx=0.01; %10 um discretization works fine
fd_param=fd2d_param(dx);


ff=3:0.5:11; % frequency of input in kHz
tau=1.28; %coefficient of the function that describes cochlear amplification, see manuscript
c=0.45;


levs=[30:20:90];
loc=360*fac; % at the 360th section in the model, CF~9 kHz
me_scale_factor=0.004; %empirical factor to scale thevenin pressure at cochlear entrance in fd2d model 
zknee=9.*(9./fd_param.fr);  %z_knee in the manuscript, which controls suppression

tic;
parfor j=1:numel(levs)
    me_g=me_gain(ff);
    vol=db2mag(levs(j))*20e-6;
    kk=0;
for fo = ff 
    kk=kk+1;
    meg=me_scale_factor*me_g(kk);
    b=fo./fd_param.fr; %beta
    dy=b*tau*(1i)-c;
    g=ones(1,numel(fd_param.x));
    for ll=1:12 %iteration for approximating nonlinear amplification in frequency domain model 
    gdy=g.*dy; %amplification
    param=param_solver(fd_param,fo,gdy);
    param.Pth=vol*meg; %thevenin pressure at cochlear entrance (i.e., model input)
    [Pd0,Vbm,Vcp,Pdbar,alpha,Pd] = fd2d_tapered(param); %change fd2d with wkb_tapered to use wkb approximation
    Vtm=2*Vcp-Vbm;
    nl=abs(Vbm./(2*pi*fo)*1e3);
    g=1./(1+(abs(nl./zknee)));
    end
    V=Vbm(loc); %vel (mm/s)
    Y=V./(2*pi*fo)*1e3; %disp (nm)?note that frequencies are in kHz
    G=Y./vol;
    sol(j).Vbm(kk)=V;
    sol(j).Ybm(kk)=Y;
    sol(j).Gbm(kk)=G;
    V=(2*Vcp(loc)-Vbm(loc)); %vel (mm/s)
    Y=V./(2*pi*fo)*1e3; %BM disp (nm)?note that frequencies are in kHz
    G=Y./vol; %gain, as reported in the experimental papers form mouse cochlea
    sol(j).Vtm(kk)=V;
    sol(j).Ytm(kk)=Y;
    sol(j).Gtm(kk)=G;
end
end
toc;
vol=db2mag(levs(end))*20e-6;
kk=0;
me_g=me_gain(ff);
tic;
for fo = ff 
    kk=kk+1;
    meg=me_scale_factor*me_g(kk);
    b=fo./fd_param.fr;
    a=0;
    param=param_solver(fd_param,fo,0); %linear model, no amplification (gdy=0)
    param.Pth=vol*meg;
    [Pd0,Vbm,Vcp,Pdbar,alpha,Pd] = fd2d_tapered(param);
    V=Vbm(loc); %vel (mm/s)
    Y=V./(2*pi*fo)*1e3; %disp (nm)?note that frequencies are in kHz
    G=Y./vol; %gain, as reported in the experimental papers form mouse cochlea
    soldead.Vbm(kk)=V;
    soldead.Ybm(kk)=Y;
    soldead.Gbm(kk)=G; %the TM solution is the same as BM post-mortem in this model
    V=(2*Vcp(loc)-Vbm(loc)); %vel (mm/s)
    Y=V./(2*pi*fo)*1e3; %BM disp (nm)?note that frequencies are in kHz
    G=Y./vol; %gain, as reported in the experimental papers form mouse cochlea
    soldead.Vtm(kk)=V;
    soldead.Ytm(kk)=Y;
    soldead.Gtm(kk)=G;

end
toc;

%% plot the transfer function
xxlim=[-1.5 0.5];
figure;
subplot(2,2,1)
[cc,id]=max(abs(sol(1).Gbm)); %find CF at the "virtual recording site" for plotting the data
cf=ff(id);
for j=1:numel(levs)
plot(log2(ff./cf),mag2db(abs(sol(j).Gbm))); hold all;
end
plot(log2(ff./cf),mag2db(abs(soldead.Gbm)),'k'); hold all;
xlim(xxlim);
xlabel('Normalized frequency (\beta)');
ylabel('BM gain [dB re 1 nm/Pa]');
set(gca,'Fontsize',12);

subplot(2,2,2);

for j=1:numel(levs)
plot(log2(ff./cf),unwrap(angle(sol(j).Ybm))/2/pi+2.); hold all;
end
plot(log2(ff./cf),unwrap(angle(soldead.Ybm))/2/pi+2.,'k');
xlim(xxlim);
xlabel('Normalized frequency (\beta)');
ylabel('BM phase [cycles]');
set(gca,'Fontsize',12);

subplot(2,2,3)
% [cc,id]=max(abs(sol(1).Gbm)); %find CF at the "virtual recording site" for plotting the data
% cf=ff(id);
for j=1:numel(levs)
plot(log2(ff./cf),mag2db(abs(sol(j).Gtm))); hold all;
end
plot(log2(ff./cf),mag2db(abs(soldead.Gbm)),'k'); hold all;
xlim(xxlim);
xlabel('Normalized frequency (\beta)');
ylabel('TM gain [dB re 1 nm/Pa]');
set(gca,'Fontsize',12);

subplot(2,2,4);

for j=1:numel(levs)
plot(log2(ff./cf),unwrap(angle(sol(j).Ytm))/2/pi+2.); hold all;
end
plot(log2(ff./cf),unwrap(angle(soldead.Ybm))/2/pi+2.,'k');
xlim(xxlim);
xlabel('Normalized frequency (\beta)');
ylabel('TM phase [cycles]');
set(gca,'Fontsize',12);
