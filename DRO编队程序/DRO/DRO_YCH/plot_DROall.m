clear
close all
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
% load('DRO_all3.mat')
load('DRO_all4.mat')
% load('DRO_all_SEsyPeriod.mat')
mu = 0.01215;
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);
T_DRO_all_day = T_DRO_all*27.2844292118811;
% for ii_loop = [1,11,21,30,40,50,60,70,80,87,91,94,96,97,98,99,100] % DRO_all2
for ii_loop = [1,2,4,6,9,12,16,24,32,48,64,96,128,192,256,300]
% for ii_loop = 1:74:300
    DT = T_DRO_all(ii_loop)*2*pi;
    t_sample = linspace(0,DT,300);
    state_ini = state_ini_all(ii_loop,:);
    sol = ode113(@(t,x)pcr3bp(t,x,mu),[0 DT], state_ini, opts_ode); 
    sol_sample = deval(sol,t_sample);
    r_DRO = sol_sample(1:2,:)'-[1-mu,0]; % 月心旋转系
%     r_DRO = sol_sample(1:2,:)'; % 共同质心旋转系
    X = [ r_DRO(:,1), r_DRO(:,2),  T_DRO_all_day(ii_loop)*ones(length(r_DRO(:,1)),1)];
    figure(1)
    patch(X(:,1),X(:,2),X(:,3),'LineWidth',1.5,'EdgeColor','interp','FaceColor','none');
    xlabel('\itx_M \rm[LU]')
    ylabel('\ity_M \rm[LU]')
    set(gca,'FontSize',15)
    axis equal; box on
end
colorbar
ylabel(colorbar,'周期 [day]','FontSize',15);
hold on
plot(0,0,'.','Color',[0.5,0.5,0.5],'MarkerSize',10)
plot(-1,0,'.','Color',[0.5,0.5,0.5],'MarkerSize',10)