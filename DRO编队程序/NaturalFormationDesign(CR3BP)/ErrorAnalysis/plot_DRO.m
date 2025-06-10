%% 惯性系和M坐标系的DRO
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../../subF_eom(CR3BP)')
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
%% 常数与变量
load('FloquetEig12')
% load('FloquetEig12_sy')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 积分
dt = 10*para.T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
sol = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 dt], x0_DRO_M_3d, opts);
abs_motion_M = deval(sol,t_sample);
abs_motion_Mclassical = -abs_motion_M; % transfer to classical M frmae
% transfer to inertial frame
abs_motion_inertial = synodic2inertial([abs_motion_Mclassical(1:3,:)+[1-con.mu;0;0]; abs_motion_Mclassical(4:6,:)],t_sample);

%% plot
figure(2)
x_ratio = 2.2; % x轴显示与图形中点的比例
y2x_ratio = 1; % 画图时y轴显示与x轴显示的比例
subplot(1,2,1); hold off
plot(abs_motion_inertial(1,:),abs_motion_inertial(2,:),'LineWidth',1.5); hold on
p1 = plot(abs_motion_inertial(1,1),abs_motion_inertial(2,1),'g^'); hold off
xlabel('\itX^{\prime} \rm[LU]'); ylabel('\itY^{\prime} \rm[LU]');
legend(p1,'初始位置')
axis equal; set(gca,'FontSize',13)
% title('DRO (inertial frame)');
title('DRO (惯性系)'); 
grid on; grid minor
x_max = max(abs_motion_inertial(1,:));
x_min = min(abs_motion_inertial(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(abs_motion_inertial(2,:))+min(abs_motion_inertial(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

subplot(1,2,2); hold off
plot(abs_motion_Mclassical(1,:),abs_motion_Mclassical(2,:),'LineWidth',1.5); hold on
p1 = plot(abs_motion_Mclassical(1,1),abs_motion_Mclassical(2,1),'g^'); hold off

xlabel('\itX \rm[LU]'); ylabel('\itY \rm[LU]');
legend(p1,'初始位置')
% legend(p1,'initial position')
axis equal; set(gca,'FontSize',13)
% title('DRO (M frame)');
title('DRO (M 坐标系)');

grid on; grid minor
x_max = max(abs_motion_Mclassical(1,:));
x_min = min(abs_motion_Mclassical(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(abs_motion_Mclassical(2,:))+min(abs_motion_Mclassical(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

% set(gcf,'Color',[255,255,255]/255);
% export_fig DRO.png -r600