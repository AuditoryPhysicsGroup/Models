%This code accompanies the manuscript by  Altoè and Charaziak (2023)
% "Intracochlear Overdrive: Characterizing Nonlinear Wave Amplification in the Mouse Apex" 
% The Journal of the Acoustical Society of America 154.5 pp 3414-3428 



%example script for calculating two-tone suppression in the model
% and generate iso-suppression tuning curves
clear all;
close all;
clc;
fac=1;
dx=0.01; %10 um discretization is a good trade-off between fine resolution and computational time
fd_param=fd2d_param(dx);

tau=1.28; %coefficient of the function that describes cochlear amplification, see manuscript
c=0.45;

loc=360; % at the 360th section in the model, CF~9 kHz
CF=9;
me_scale_factor=0.002; %empirical factor to scale thevenin pressure at cochlear entrance in fd2d model 
ss=9.*(9./fd_param.fr); %z_knee in the manuscript, which controls suppression

probe_frequency=9.1;
supp_frequency=1:0.5:18;
levs=[30:5:90]; % suppressor_levels;
ff=1:0.5:15; % suppressor frequency

probe_level=40; % probe level in dB probe;
vol_probe=db2mag(probe_level)*20e-6;
Pprobe=me_scale_factor*me_gain(probe_frequency)*vol_probe;
bprobe=probe_frequency./fd_param.fr;
dy_probe=bprobe*tau*(1i)-c;
g=ones(1,numel(fd_param.x));
for ll=1:12 %iteration for approximating nonlinear amplification in frequency domain model 
    param=param_solver(fd_param,probe_frequency,g.*dy_probe);
    param.Pth=Pprobe; %thevenin pressure at cochlear entrance (i.e., model input)
    [Pd0,Vbm,Vcp,Pdbar,alpha,Pd] = wkb_tapered(param);
    nl=abs(Vbm./(2*pi*probe_frequency)*1e3);
    g=1./(1+(abs(nl./ss)));
end
Vprobe=Vbm(loc);
Yprobe_alone=Vprobe/(2*pi*probe_frequency)*1e3;

clear Vbm Vcp Pdbar alpha Pd
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
    param=param_solver(fd_param,fo,g.*dy); % responses at suppressor frequency
    param.Pth=vol*meg; %thevenin pressure at cochlear entrance (i.e., model input)
    [Pd0,Vbm_sup,Vcp_sup,Pdbar,alpha,Pd] = wkb_tapered(param);
    param=param_solver(fd_param,probe_frequency,g.*dy_probe); % responses at probe frequency
    param.Pth=Pprobe; %thevenin pressure at cochlear entrance (i.e., model input)
    [Pd0,Vbm_probe,Vcp_probe,Pdbar,alpha,Pd] = wkb_tapered(param);
    %overall nonlinearity rms BM displacement
    nl=sqrt(abs(Vbm_probe./(2*pi*probe_frequency)*1e3).^2+abs(Vbm_sup./(2*pi*fo)*1e3).^2);
    g=1./(1+((nl./ss)));
    end
    Vp=Vbm_probe(loc); %vel (mm/s)
    Yp=Vp./(2*pi*probe_frequency)*1e3; %disp (nm)?note that frequencies are in kHz
    Vs=Vbm_sup(loc); %vel (mm/s)
    Ys=Vs./(2*pi*fo)*1e3; %disp (nm) note that frequencies are in kHz
    sol(j).Ys(kk)=Ys;
    sol(j).Yp(kk)=Yp;
end
end
toc;

for j=1:numel(levs)
   Yprobe(j,:)=sol(j).Yp;
   Ysup(j,:)=sol(j).Ys;
end
%% plot iso suppression tuning curves
supp_db=mag2db(abs(Yprobe_alone./Yprobe));
threshold=6; %suppression threshold in dB

   for i=1:numel(ff)
       tmp=supp_db(:,i)-threshold; %subtract threshold from suppression in dB
       prod=tmp(1:end-1).*tmp(2:end); %find the two levels across which the data meets the threshold critetion 
       idxs=find(sign(prod)<0);
       if(numel(idxs))
       id=idxs(1)
       iso_sup(i)=interp1(tmp(id-1:id+1),levs(id-1:id+1),0,'linear'); % linera interpolate to find the suppression level (as in Dewey 2019)
       sup_disp(i)=db2mag(interp1(levs(id-1:id+1),mag2db(abs(Ysup((id-1:id+1),i))),iso_sup(i),'pchip'));
       else
       iso_sup(i)=nan;
        sup_disp(i)=nan;
       end
   end

figure(1);
plot(log2(ff/(CF)),iso_sup); hold all;
xlabel('Suppressor frequency [octave re CF]');
ylabel('Suppressor level [dB SPL]');
xlim([-2.5 0.5]);


