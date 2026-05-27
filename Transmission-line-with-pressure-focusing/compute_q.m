function q=compute_q(model,t)
fs=model.fs;
tsample=t*fs;
ta=floor(t*fs)+1;
a=model.stim(ta);
b=model.stim(ta+1);
c=model.stim(ta+2);
d=model.stim(ta+3);

frac=tsample-floor(tsample);
cminusb=c-b;
stim=b + frac .* ...
        (cminusb - 0.1666667 .* (1. - frac) .*...
         ((d - a - 3.0 .* cminusb) .* frac + (d + 2.0 .* a - 3.0 .* b)));
model.q(1)=stim./model.mme;
q=model.q;
