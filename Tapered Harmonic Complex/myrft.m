function [ft,fax]=myrft(sig,fs)

    ft=rfft(sig);
    fax=1:1:numel(ft);
    fax=(fax-1);
    fax=fax/fax(end)*fs/2;