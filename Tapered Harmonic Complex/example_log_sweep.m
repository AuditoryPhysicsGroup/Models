% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)
% this code plots Figure 1 of that manuscript

fs=96e3;
r=20; %rep rates


f1=1.6e3;%1.6e3; %low cut-off
f2=6.4e3; %high cut-off

[freq,rate,T]=sweepPreliminaries(r,fs); %adjust rate and frequencies so signals get periodic in digital domain
dur=1/rate;
Nsamples=round(dur*fs);
Nrep=20;

order=16;
AS=log_tapered_complex_amp(freq,f1,f2,order); % define amplitude
[specS,phaseS]=make_complex_tone_spectrum(AS); %use schroeder's equation to generate spectrum
toneS=irfft(specS);
toneS=toneS./max(toneS);
soundS=repmat(toneS,1,Nrep);
t=((1:1:numel(soundS))-1)/fs;


figure; subplot(221); plot(t*1e3,soundS); xlabel('Time [ms]'); ylabel('Amplitude');
xlim([0 200]);
subplot(222); plot(freq,mag2db(abs(specS))); xlabel('Frequency [Hz]'); ylabel('Amplitude [dB]');
xlim([500 10e3]); ylim([-80 5]);
subplot(2,2,3:4); spectrogram(soundS(1:4*fs/rate),256,255,512,fs,'yaxis');
ylim([0.1 10]);
  clim([-85 -35]);