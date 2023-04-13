%
% Mechanical Toy Model of OoC 
%
% Accompanying Altoè A, Dewey JB, Charaziak KK, Oghalai JS, Shera CA. 
%    Overturning the mechanisms of cochlear amplification via area deformations of the organ of Corti.
%    J. Acoust. Soc. Am. 152 2022
% See Appendix A
%

close all;
f = 0.01:0.01:2;                        % frequency axis
s = 1i*2*pi*f;                          % Laplace s (=i*w=i*2pi*f)

Kdc = 1;                                % stiffness of Deiters Cell
Ktm = 2*Kdc;                            % stiffness of RL/TM complex
Kbm = 10*Kdc;                           % stiffness of BM
C = 1;                                  % damping coefficient of Deiters cells
Zbm = Kbm;                              % BM impedance just stiffness
Ztm = Ktm;                              % RL/TM impedace just stiffness

% Deiters Cell impedance parallel of spring and damper
Zdc = Kdc + s*C;

% equivalent combined impedance of BM and DC (i.e., mechanical load at
% bottom of OHC) spring in series with Zdc.  [spring connected in series
% with DC spring-damper impedance.]
Zeq = Zbm.*Zdc./(Zbm + Zdc); 

% motion of the OHC/DC relative to BM (following unitary
% contraction/expansion of OHC)
tf_ohc_bottom = Ztm./(Zeq + Ztm); 

% motion of the apical surface of the OHC
tf_ohc_top = -Zeq.*tf_ohc_bottom./Ztm; 

% motion of the BM.
tf_bm = -tf_ohc_top.*Ztm./Zbm; 

figure; 
subplot(1,2,1);
semilogx(f,mag2db(abs(tf_ohc_top)));
hold all; 
semilogx(f,mag2db(abs(tf_ohc_bottom)));
semilogx(f,mag2db(abs(tf_bm)));
legend('TM','DC','BM');
xlabel('Normalized Frequency');
ylabel('Displacement re OHC [dB]');
set(gca,'fontsize',11);

% test: if calculations the following plot should be a straight line
% with value 1 (OHC expansion = OHC_top-OHC_bottom = 1).
% figure; semilogx(f,(abs(tf_ohc_top-tf_ohc_bottom))); 

subplot(1,2,2);
semilogx(f,wrapToPi(angle(tf_ohc_top))/(2*pi));
hold all; 
semilogx(f,wrapToPi(angle(tf_ohc_bottom))/(2*pi));
semilogx(f,wrapTo2Pi(angle(tf_bm))/(2*pi));
legend('TM','DC','BM');
xlabel('Normalized Frequency');
ylabel('Phase re OHC [cyc]');
set(gca,'fontsize',11);
% print('-depsc2','toy model')

%
% Superimpose overall translation of the OOC and plot mag and phase
%
translation = -1i*tf_ohc_top/2;
% motion of the OHC/DC relative to BM (following unitary
% contraction/expansion of OHC)
tf_ohc_bottom_s = tf_ohc_bottom + translation; 
% motion of the apical surface of the OHC
tf_ohc_top_s = tf_ohc_top + translation; 
tf_bm_s = tf_bm + translation;

Ttop_bm = tf_ohc_top_s./tf_bm_s;
Tbot_bm = tf_ohc_bottom_s./tf_bm_s;

figure; 
subplot(1,2,1);
semilogx(f,mag2db(abs(Ttop_bm)));
hold all; 
semilogx(f,mag2db(abs(Tbot_bm)));
legend('TM','DC');
xlabel('Normalized Frequency');
ylabel('Displacement re BM [dB]');
set(gca,'fontsize',11);
% print('-depsc2','Mechanical_model');

subplot(1,2,2);
semilogx(f,unwrap(angle(Ttop_bm))/(2*pi));
hold all; 
semilogx(f,unwrap(angle(Tbot_bm))/(2*pi));
xlabel('Normalized Frequency');
ylabel('Phase re BM [cyc]');
set(gca,'fontsize',11);
% print('-depsc2','relative_motions')
