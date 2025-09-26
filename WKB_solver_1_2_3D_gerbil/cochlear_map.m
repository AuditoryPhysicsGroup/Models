function [x,fr,mapParams] = cochlear_map (dx,species,flag)
% Computes cochlear map or its inverse for several species.
% function [fr,mapParams] = cochlear_map (x_over_L,?species='human'?,?flag=1?)
%
% If flag>0, returns cochlear position-frequency map fr(x_over_L)
% in Hz, where x_over_L is normalized distance from the base
% (i.e, x/L, where L is the length the cochlea and 0<=x/L<=1).
%
% If flag<0, inverts the operation and returns x_over_L(fr):
% function x_over_L = cochlear_map (fr,?species='human'?,-1)
%
% Based on Greenwood (1990) curves for several species.  Supported
% species are 'human', 'cat', 'gerbil', 'guinea pig', 'guinea pig
% greenwood', 'chinchilla', 'chinchilla muller', 'chinchilla muller
% shifted', 'mouse', and 'erb' (meaning based on data from Glasberg
% and Moore 1990 assuming an ERB corresponds to a constant distance
% along the BM). See also greenwood() and cochmap(). Default guinea
% pig data now from Tsuji and Liberman.
  
% Conversion from Greenwood form to (f0+f1)e^(-x/d)-f1
%   f0 = A*(10^a-k)
%   f1 = A*k
%   d = L/(a*log(10))

if (nargin < 2 | isempty(species)), species = 'human'; end
if (nargin < 3 | isempty(flag)), flag = 1; end

switch (species)
 case {'human'}
  A = 165.4;
  L = 35;
  a = 2.1;
  k = 1;
 case {'human-exponential'}
  A = 165.4;
  L = 35;
  a = 2.1;
  k = 0;
 case {'erb'}
  A = 165.4;
  L = 35;
  a = 2.1;
  k = A*4.37/1000;			% see Glasberg and Moore
                                        % 1990, Hearing Res., 47:103-138
 case {'cat'}
  A = 456;
  L = 25;
  a = 2.1;
  k = 0.8;
 case {'guinea pig greenwood'}
   A = 350;
   L = 18.5;
   a = 2.1;
   k = 0.85;
 case {'guinea pig'}
  % From Tsuji and Liberman 1997 J. Comp. Neurol. 381:188-202
  A = 131.9;
  L = 20;
  a = 2.618;
  k = 0.0;
 case {'chinchilla'}			% from Greenwood 1990
   A = 163.5;
   L = 18.4;
   a = 2.1;
   k = 0.85;
 case {'chinchilla muller'}             % exponential from Muller et al HR 2010
   A = 120;
   L = 20;
   a = 2.37;
   k = 0.0;
 case {'chinchilla muller shifted'}     % shifted exponential from Muller et al HR 2010
   A = 98.9;
   L = 20;
   a = 2.47;
   k = 0.685;
 case {'mouse'}				% from Muller et al HR 2005
   A = 4488;
   L = 5.13;
   a = 1.2;
   k = 0.0;
 case {'gerbil'}			% from Muller HR 1996
   A = 398;
   L = 11.1;
   a = 2.2;
   k = 0.631;
 otherwise
  disp ('unsupported species');
  fr = [];
  return
end

  x=0:dx:L;
  x_over_L=x/L;
  fr = A*(10.^(a*(1-x_over_L))-k);
  f0 = A*(10^a-k);
  f1 = A*k;
  d = L/(a*log(10));

  mapParams.species = species;
  mapParams.A = A;
  mapParams.L = L;
  mapParams.a = a;
  mapParams.k = k;
  mapParams.fmax = f0;
  mapParams.fcut = f1;
  mapParams.d = d;
end




