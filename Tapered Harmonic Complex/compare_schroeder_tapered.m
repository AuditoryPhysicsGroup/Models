% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)
% this code plots Figure 1 of that manuscript


fs=96e3;
r=50; %rep rate


f1=1.6e3;%1.6e3; %low cut-off
f2=6.4e3; %high cut-off

[freq,rate,T]=sweepPreliminaries(r,fs); %adjust rate and frequencies so signals get periodic in digital domain
dur=1/rate;
Nsamples=round(dur*fs);

%periodic FM sweep
[toneFM]=periodic_FM_sweep(f1,f2,rate,fs);
toneFM=toneFM*sqrt((f2-f1)*T);

% Schroeder's harmonic complex
AS=schroeder_amp(freq,f1,f2); % define amplitude
[specS,phaseS]=make_complex_tone_spectrum(AS); %use schroeder's equation to generate spectrum
toneS=irfft(specS); %generate tone via ifft

% Tapered Harmonic Complex
order=16; % order of filtering per tapering
AT=linear_tapered_complex_amp(freq,f1,f2,order); % define amplitude
[specT,phaseT]=make_complex_tone_spectrum(AT); %use schroeder's equation to generate spectrum from amplitude
toneT=irfft(specT); %generate tone via ifft
t=(0:1:(numel(toneT)-1))/fs;


%% repeat and normalize
Nrep=20;
soundSch=repmat(toneS,1,Nrep);

soundT=repmat(toneT,1,Nrep);

soundFM=repmat(toneFM,1,Nrep);

% 
soundFM=soundFM./max(soundFM); 
soundT=soundT./max(soundT);
soundSch=soundSch/max(soundSch);
%uncomment to save files
% audiowrite('Tapered_Schoreder_1.6_6.4_kHz_50Hz.wav',soundT*0.9,fs);
% audiowrite('Regular_Schoreder_1.6_6.4_kHz_50Hz.wav',soundSch*0.9,fs);
% audiowrite('FMSweep_1.6_6.4_kHz_50Hz.wav',soundFM*0.9,fs);


%%
figure; subplot(3,3,3);
spectrogram(soundSch(1:4*fs/rate),256,255,512,fs,'yaxis');

% title('Schroeder complex');

subplot(3,3,6);
spectrogram(soundFM(1:4*fs/rate),256,255,512,fs,'yaxis');
% title('FM sweep')


subplot(3,3,9);
spectrogram(soundT(1:4*fs/rate),256,255,512,fs,'yaxis');
% title('Tapered sweep')


tms=((1:1:fs/rate)-1)/fs*1e3

 subplot(3,3,1);
 tms=((1:1:numel(soundSch))-1)/fs*1e3;
plot(tms,soundSch,'k','linewidth',0.75);
% clim([-40 0]);
% ylim([0.05 10]);
% title('Schroeder complex');

subplot(3,3,4);
 tms=((1:1:numel(soundFM))-1)/fs*1e3;

plot(tms,soundFM,'k','linewidth',0.75);



subplot(3,3,7);
tms=((1:1:numel(soundT))-1)/fs*1e3;
plot(tms,soundT,'k','linewidth',0.75);



subplot(3,3,2);
plot(freq/1e3,AS,'k','linewidth',2);
xlim([0 10]);
ylim([0 1.2]);
title('Schroeder Complex');


subplot(3,3,8);
plot(freq/1e3,AT,'k','linewidth',2);
xlim([0 10]);
ylim([0 1.2]);
title('Tapered Sweep')


subplot(3,3,5);
[ft,fx]=myrft(soundFM(1:Nsamples),fs);
plot(fx/1e3,abs(ft)./max(abs(ft)),'k','linewidth',2);
xlim([0 10]);
ylim([0 1.2]);
title('FM Sweep')

%%

for i=[2,5,8]
    subplot(3,3,i);
xlims=[0 10];
xlabel('Frequency [kHz]');
ylabel('Amplitude');
ylim([0 1.2]);
set(gca,'fontsize',12);
xlim(xlims);
box 'off'
end
plot_rep=2/rate*1e3;
xlims=[0 plot_rep];
for i=[1,4,7]
    subplot(3,3,i);
xlim(xlims);
ylim([-1.2 1.2])
box 'off';
xlabel('Time [ms]');
ylabel('Amplitude');
set(gca,'FontSize',12);
end

plot_rep=4/rate*1e3;
xlims=[0 plot_rep];
ylims=[0 12];

climss=[-85 -35];
for i=[3,6,9]
    subplot(3,3,i);
    clim(climss);
    ylim(ylims);
    xlim(xlims);
    set(gca,'FontSize',12);
    box 'off'

end

