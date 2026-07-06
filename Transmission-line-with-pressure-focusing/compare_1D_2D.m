
clear all;
close all;
clc;

figure(77); 

fs=100e3;

dx=0.02;
f0=16e3;

iteration=0;
for tau=[0 0.55] % we start with a passive model and then repeat calculation for active model
[model1D,model2D]=gerbil_state_space(dx,tau);
param=gerbil_box_model(dx);

t=0:1/fs:0.05; % 50 ms sound
tt=0.25e-3;%-3;
stim=ramp(sin(2*pi*f0*t),tt*fs);

model=model2D;
model.t=t;
model.stim=set_model_stim(stim.');%%
model.fs=fs;
tic;
% compute the time domain solution
options = odeset('RelTol',1e-3,'AbsTol',1e-9);
[~,yy]=ode45(@(t,y) state_space_solver_fast2D(t,y,model),t,(model.ystart),options);
toc;
vbm=yy(:,3:4:end-1); %extract the solution of vbm and vint 
vint=yy(:,5:4:end-1);
for i=1:model.Nx
    Vt2d(i)=fitKnownSines(t,-vbm(:,i),f0); % calculate amplitude and phase of BM response at the stimulus frequency
    % the minus sign is because in the model positive velocity towards
    % scala tympani (experimenter typically show the opposite)
end
vbm2D=vbm;

tic;
model=model1D;
model.t=t;
model.stim=set_model_stim(stim.');%%
model.fs=fs;
tic;
% options = odeset('RelTol',1e-3,'AbsTol',1e-9);
[~,yy]=ode45(@(t,y) state_space_solver_fast1D(t,y,model),t,(model.ystart),options);
toc;
vbm=yy(:,3:4:end-1);
vint=yy(:,5:4:end-1);
for i=1:model.Nx
    Vt1d(i)=fitKnownSines(t,-vbm(:,i),f0);
end

% calculate passive Finite difference solution
params=param_solver(param,f0,tau);
[~,VbmFE]=fd2d_tapered(params);


figure(77);
   subplot(221);
if (iteration==0)
 
plot(param.x,mag2db(abs(Vt2d./Vt2d(1))),'r:'); hold all; 
plot(param.x,mag2db(abs(VbmFE./VbmFE(1))),'k:'); 
plot(param.x,mag2db(abs(Vt1d./Vt1d(1))),'b:'); 

%normalize the curves, the input pressure is defined slightly different in the finite difference and time domain solver  
else
set(gca,'FontSize',14);
% xlabel('Distance from Cochlear entrance [mm]');
% ylabel('BM velocity [dB]');
ylim([-30 50])
box off;

plot(param.x,mag2db(abs(Vt2d./Vt2d(1))),'r'); hold all; 
plot(param.x,mag2db(abs(VbmFE./VbmFE(1))),'k'); 
plot(param.x,mag2db(abs(Vt1d./Vt1d(1))),'b'); 

end
% plot(param.x,mag2db(abs(Vt1d./ctd)),'r:'); 
% plot(param.x,mag2db(abs(Vt2d./ctd)),'b:'); 

xlim([0 4]);
xlabel('Distance along the BM [mm]')
ylabel('BM velocity magnitude [dB]')


 % hold(ax(2),'on');  
subplot(222)
if (iteration==0)


plot(param.x,unwrap(angle(VbmFE./VbmFE(1)))/2/pi,'k:'); hold all; 
plot(param.x,unwrap(angle(Vt2d./Vt2d(1)))/2/pi,'r:');
plot(param.x,unwrap(angle(Vt1d./Vt1d(1)))/2/pi,'b:');

else

plot(param.x,unwrap(angle(VbmFE./VbmFE(1)))/2/pi,'k'); hold all; 
plot(param.x,unwrap(angle(Vt2d./Vt2d(1)))/2/pi,'r');
plot(param.x,unwrap(angle(Vt1d./Vt1d(1)))/2/pi,'b');
end 
xlim([0 4]);
ylim([-3 0])
set(gca,'FontSize',14);
xlabel('Distance along the BM [mm]')
ylabel('BM velocity phase [cycle]')
iteration=iteration+1;
end

