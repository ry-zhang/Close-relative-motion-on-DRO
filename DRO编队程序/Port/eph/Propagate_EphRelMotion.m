clear
addpath('./subF_eom(eph)')

format longg
format compact
warning off

%% ga优化星历下的自然编队构型
aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC  = [2030 1 1 0 0 0]; % 初始历元
aux = initialize(aux); % 初始化

load('xTC_VVLH_FHL_59days.mat')
aux.jd0  = jd0; % 初始历元
aux.t0UTC = [];

% 星历积分区间 
tspan_sec = [0,29.5*86400*2];% 
t_sample = linspace(tspan_sec(1),tspan_sec(2),500);
% 主星星历积分
[xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);

%% 计算优化结果
% 副星初值转换至地月旋转系
x0_MCR_chaser = T_TCO2TCR_eph(xTC_VVLH,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);
% 计算相对运动
rvTC_MCR = xx_MCR_chaser-xx_MCR_target;
rvTC_VVLH = T_TCR2TCO_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'VVLH');

% 画相对运动
figure(1)
plot(rvTC_VVLH(:,1), rvTC_VVLH(:,3),'color',[0 0.4470 0.7410],'LineWidth',2); hold on;
plot(rvTC_VVLH(1,1), rvTC_VVLH(1,3),'g^');
plot(rvTC_VVLH(end,1), rvTC_VVLH(end,3),'rv');
box on; grid on; grid minor; 
hold off; 
axis equal; xlabel('x[km]'); ylabel('z[km]')
set(gca,'FontSize',15,'fontname','times new roman'); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title('RelMotion(VVLH)')

%% 画DRO
f1 = figure(2);
set(f1,'name','星历DRO MCR')
% figure('color',[1 1 1],'name','星历DRO MCR')
p2 = plot(xx_MCR_target(:,1),xx_MCR_target(:,2),'color',[0 0.4470 0.7410],'LineWidth',1.5); hold on;
plot(xx_MCR_target(1,1),xx_MCR_target(1,2),'g^');
plot(xx_MCR_target(end,1),xx_MCR_target(end,2),'rv');
box on; grid on; grid minor; hold off;
axis equal; xlabel('x[km]'); ylabel('y[km]')
set(gca,'FontSize',15,'fontname','times new roman');
title('DRO')
