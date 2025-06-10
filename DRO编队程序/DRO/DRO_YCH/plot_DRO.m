clear
close all
set(0,'defaultAxesFontName', 'TimesSimSun','defaultTextFontName', 'TimesSimSun');
set(0,'defaultAxesFontSize',15,'defaultTextFontSize',15)
set(0,'defaultLineLineWidth',1.5)

%% 常数与变量
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2
% con.mu = 0.01211; % 
con.mu = 0.01215; % 20200531
con.r_norma = 3.84399*10^5; % km
% T_M = 27.321661; % day
con.T_norma = sqrt(con.r_norma^3/(mu_E+mu_M)); % s
con.T_norma_day = con.T_norma/3600/24;
con.v_norma = con.r_norma/con.T_norma; % km/s
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 画图
% load('DRO_all3.mat')
% load('DRO_all_SEsyPeriod.mat')
load('DRO_all_MoonPeriod.mat')
% T_DRO_all = dataDRO(:,7)/pi;
% state_ini_all = dataDRO(:,[1,2,4,5]);
% for ii_loop = 1:2:length(T_DRO_all)
f1 = figure(1);
set(f1,'Units','centimeters','Position',[35,15,10,8])
f2 = figure(2);
set(f2,'Units','centimeters','Position',[35,5,10,8])
for ii_loop = [4,6,7]
% for ii_loop = [4]
    DT = T_DRO_all(ii_loop)*2*pi;
    t_sample = linspace(0,4*DT,500);
    state_ini = state_ini_all(ii_loop,:);
    sol = ode113(@(t,x)pcr3bp(t,x,con.mu),[0 4*DT], state_ini, opts_ode); 
    sol_sample = deval(sol,t_sample);
    r_DRO_M = sol_sample(1:2,:)'-[1-con.mu,0]; % 月心旋转系
    r_DRO_I = synodic2inertial(sol_sample,t_sample)';
    
    figure(1)
    pO = plot(r_DRO_M(:,1)*con.r_norma,r_DRO_M(:,2)*con.r_norma,'LineWidth',1.5); hold on
    xlabel('\itx_M \rm[km]')
    ylabel('\ity_M \rm[km]')
%     title('DRO (月心坐标系M)')
    set(gca,'FontSize',13)
    axis equal; box on
    grid on
    
    figure(2)
    plot(r_DRO_I(:,1)*con.r_norma,r_DRO_I(:,2)*con.r_norma,'LineWidth',1.5); hold on
    xlabel('\itx_{I} \rm[km]')
    ylabel('\ity_{I} \rm[km]')
%     title('DRO (月心坐标系M)')
    set(gca,'FontSize',13)
    axis equal; box on
    grid on
end
hold on
figure(1)
pE = plot(-con.r_norma, 0,'.','Color',[5, 102, 176]/255,'MarkerSize',20);
pM = plot(0, 0,'.','Color',[133, 133, 133]/255,'MarkerSize',10);
% legend('2:1','3:1','4:1','Location','north')
hold off
% ylim([-11e4,11e4]); legend([pM,pO],{'月球','2:1 DRO'},'Location','northeastoutside');
set(gcf,'Color',[255,255,255]/255);
% export_fig DROM.png -r600

figure(2)
plot(sin(t_sample)*con.r_norma,cos(t_sample)*con.r_norma,'Color',[133, 133, 133]/255,'MarkerSize',10);
plot(0, 0,'.','Color',[5, 102, 176]/255,'MarkerSize',20);
% legend('2:1','3:1','4:1','月球轨道','Location','north')
ylim([-5e5,5e5])
hold off
set(gcf,'Color',[255,255,255]/255);
% export_fig DROI.png -r600

