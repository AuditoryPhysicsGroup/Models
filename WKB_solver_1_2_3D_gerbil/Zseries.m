function z = Zseries(param,f)
% function z = Zseries
% The values of unspecified arguments are
% taken from corresponding global variables

% S=param.H.^2*pi/2;
%  S0=pi*param.H(1)^2/2;
%  S=S0.*param.acoustic_area_tapering;
 z = 1i*2*pi*f*param.M;

% z = 1i*beta*wmaxL1;			% tapered
% if (~tapered)
%   z = z.*(fr/fmax);                  % same as i*w*L1
% end


