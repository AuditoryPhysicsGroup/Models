close all; clear all; clc;

%    Example file on how to use the inner-hair-cell model "IHC_model.m", supplementing 
%    Estimation of the effect of the voltage activation of the different
%    basolateral K+ channels on the IHC DC response to a sinusoidal
%    stimulus

%    Script supplementing:
%    "Dierich et al., Optimized Tuning of Auditory Inner Hair Cells to Encode Complex Sound through 
%    Synergistic Activity of Six Independent K+ Current Entities", Cell Reports (2020) 32, 107869
%    https://doi.org/10.1016/j.celrep.2020.107869
%    For fair use only.

fs=50e3; % sample frequency
f0=1e3; % frequency of input to model
t=0:1/fs:55e-3; %time axis
ton=1e-3*fs; %onset/offset ramp
tend1=50e-3*fs; %length of sinusoidal stimulation

sig=sin(2*pi*f0*t); %input sinusoid
sig(1:ton)=sig(1:ton).*linspace(0,1,ton);
sig(tend1+1:end)=0;
sig(tend1-ton+1:tend1)=sig(tend1-ton+1:tend1).*linspace(1,0,ton);
amplitude=logspace(log10(5),log10(500),30)*1e-9;
for jj=1:numel(amplitude)
    
ihcIn=sig*amplitude(jj); %input of the model is deflection of stereocilia bundle (specified in meter)

Vnative=IHC_model(ihcIn,fs); %use IHC_model(input,sample frequency) to run the native model (all K+ channels are voltage dependent)

% Tha last optional parameter specify which channel(s) should be made voltage-independent (constant conductance)
VnoF=IHC_model(ihcIn,fs,'Kf'); %render voltage-independent the fast-activating K+ channels

VnoS=IHC_model(ihcIn,fs,{'K11','K12','K18','K74'}); %render voltage-independent the slow-activating K+ channels


Vind=IHC_model(ihcIn,fs,{'K11','K12','K18','K74','Kf'}); %model with all voltage-independent K+ currents, as in classic IHC models


% compute DC response (DC=mean value of receptor potential)
DC_no_F(jj)=mean(VnoF-VnoF(1))*1e3; % the membrane potential at initial instant is the resting potential
DC_no_S(jj)=mean(VnoS-VnoS(1))*1e3; % note that it is identical in all simulation
DC_native(jj)=mean(Vnative-Vnative(1))*1e3;
DC_ind(jj)=mean(Vind-Vind(1))*1e3;

end

% plot the responses
figure;
semilogx(amplitude*1e9,DC_native,'k');
hold all;
semilogx(amplitude*1e9,DC_no_S);
hold all;
semilogx(amplitude*1e9,DC_no_F);
semilogx(amplitude*1e9,DC_ind);
xlim([5,500]);
xlabel('Vibration amplitude [nm]');
ylabel('DC response [mV]');
leg=legend('Native','I_{K,f} only','I_{K,s} only','Voltage independent');
set(leg,'FontSize',14,'Location','Northwest');