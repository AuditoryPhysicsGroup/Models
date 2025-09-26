clear all;
close all;
clc;
dx=0.01;
param=parameters_gerbil(dx);

flist=[1 2 4 8 16 32]*1e3; %stimulus frequencies (in Hz)
tau=1; %tau=1 active model
%tau =0; %tau=0 passive model

for i=1:numel(flist)
    ff=flist(i);
     b=ff./param.fr; %beta
    dy=b*tau*(1i); %you can try more fantastic forms for tau
    param_s=param_solver(param,ff,dy);
    %there are several wkb solvers to choose from:
    % 1D, 2D_reduced, 2D_augmented, 3D_cylindrical, 3D_box

    [k,Vbm,~,~,~]=wkb_3D_cylindrical(param_s); %3D "cylindrical" model
    sol(i).Vbm=Vbm;
 
end

cc=max(abs(sol(1).Vbm)); %normalization factor to choose a decent y-axis range 

for i=1:numel(flist)
plot(param.x,mag2db(abs(sol(i).Vbm)./cc),'k','LineWidth',1); hold all;
end
ylim([-30 40]);
set(gca,'fontsize',14);
xlabel('distance from stapes [mm]');
ylabel('BM normalized velocity [dB]')