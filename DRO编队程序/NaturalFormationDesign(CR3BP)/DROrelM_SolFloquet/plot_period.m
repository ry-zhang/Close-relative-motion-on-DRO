%% 载入数据
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
addpath('../../subF_eom(CR3BP)')
load('FloquetEig_all4.mat')
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 画图
m2n_all = 2*pi./alpha_all;
T0_all = J_period_all(:,2)*con.T_norma_day*2*pi; % day
T12_all = m2n_all.*T0_all;
omega0_all = 2*pi./T0_all;
omega12_all = 2*pi./T12_all;

figure(1)
subplot(1,2,1)
plot(T0_all,T12_all(:,1),'LineWidth',1.5); 
hold on; plot(T0_all,T12_all(:,2),'LineWidth',1.5); hold off
set(gca,'FontSize',13)
grid on; grid minor
xlabel('\itT\rm_{0} [day]')
% ylabel('\itT\rm_{1} [day]')
ylabel('周期 [day]')
legend('\itT\rm_{1}','\itT\rm_{2}','Location','north')

subplot(1,2,2)
plot(T0_all,T12_all(:,1)./T0_all,'LineWidth',1.5); hold on;
plot(T0_all,T12_all(:,2)./T0_all,'LineWidth',1.5); 
% plot([0;T0_all(1)],[1;1]*3,'--','Color',0.5*[1,1,1],'LineWidth',1.5); 
% plot([0;T0_all(1)],[1;1]*4,'--','Color',0.5*[1,1,1],'LineWidth',1.5); 
% plot([0;T0_all(1)],[1;1]*5,'--','Color',0.5*[1,1,1],'LineWidth',1.5); 
% % plot([5;T0_all],[1;ones(size(T0_all))]*6,'--','Color',0.5*[1,1,1],'LineWidth',1.5); 
hold off

set(gca,'FontSize',13)
grid on; grid minor
xlabel('\itT\rm_{0} [day]')
% ylabel('\itT\rm_{1}/\itT\rm_{0}')
ylabel('周期比')
legend('\itT\rm_{1}/\itT\rm_{0}','\itT\rm_{2}/\itT\rm_{0}','Location','north')

% set(gcf,'Color',[255,255,255]/255);
% export_fig DRO_T12.png -r600

% figure(2)
% plot3(T0_all,T12_all(:,1),T12_all(:,2),'LineWidth',1.5);
% % plot3(omega0_all,omega12_all(:,1),omega12_all(:,2),'LineWidth',1.5);
% % axis equal
% grid on; box on
% xlabel('\itT\rm_{\omega_0} [day]')
% ylabel('\itT\rm_{\omega_1} [day]')
% zlabel('\itT\rm_{\omega_2} [day]')