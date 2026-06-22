function A=linear_tapered_complex_amp(freq,f1,f2,order)
% This script accompanies the manuscript
% Altoè A., "Suppressing spectral edge effectsin Schroeder's Harmonic Complex" (2026)

% Calculate the amplitude spectrum of the Tapered Schroeder Complex


A=zeros(size(freq));
twoN=order*2;
norm=2*(1-1/(twoN^2-1)); %small correction to keep the ERB band the same as that of a Schroeder's complex
deltaf=(f1/(twoN+1)+f2/(twoN-1))/norm; %the formula is obtained by assuming the order is large
f1=f1+deltaf;
f2=f2-deltaf;
for j=2:numel(freq)-1
    ff=freq(j);
    if ff<f1
        A(j)=(ff/f1)^order; 
    elseif ff>f2
        A(j)=(f2/ff)^order;
    else
        A(j)=1;
    end
end

