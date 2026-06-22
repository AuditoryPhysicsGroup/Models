% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)
% this code plots Figure 1 of that manuscript


% create a Tapeded Schroeder Phase Complex Lentz and Leek (2001) style.
% the paramter "C" [-1,1] controls phase curvature and crest factor;
fs=96e3;
r=50; %rep rate


f1=1.6e3;%1.6e3; %low cut-off
f2=6.4e3; %high cut-off

[freq,rate,T]=sweepPreliminaries(r,fs); %adjust rate and frequencies so signals get periodic in digital domain
dur=1/rate;
Nsamples=round(dur*fs);

   order=16; % order of filtering per tapering
    AT=linear_tapered_complex_amp(freq,f1,f2,order); % define amplitude
    BT=schroeder_amp(freq,f1,f2);
C=[-1,-0.25,0.25,1];
for i=1:numel(C)  
    [specT,phaseT]=phase_complex_tone_spectrum(AT,C(i)); %use schroeder's equation to generate spectrum from amplitude
    tapered_sch=irfft(specT); %generate tone via ifft
    tapered_sch=tapered_sch./max(tapered_sch);
    subplot(2,2,i);
    plot(tapered_sch); hold all;
    [specT,phaseT]=phase_complex_tone_spectrum(BT,C(i)); 
     sch=irfft(specT); %generate tone via ifft
    sch=sch./max(sch);
        plot(sch);
    if(i==1)
        legend('Tapered Schroeder','Classic Schroeder');
    end
    title(['C=' num2str(C(i))]);
end