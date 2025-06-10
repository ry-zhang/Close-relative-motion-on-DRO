clear

set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultTextFontSize', 14)
set(0, 'DefaultLineLineWidth', 1)
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
%% 星历积分
addpath('../subF_eom(CR3BP)')
addpath('../subF_eom(eph)')

aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC  = datetime([2030 1 1 0 0 0]); % 初始历元
aux = initialize(aux); % 初始化

% 星历积分区间 
tspan_sec = [0,29.5*86400*0.5];% s
% 星历采样点
t_sample = linspace(tspan_sec(1),tspan_sec(2),100000);

% 星历积分初值（旋转系）
% x0_MCR_target = [72687.2175909459 0 1e4 0 -0.540957166578942 0]; % 星历DRO初值
x0_MCR_target = [72687.2175909459 0 0 0 -0.540957166578942 0]; % 星历DRO初值
% x0_MCR_target = [0.162164232531584*aux.LU,0,0.0791678086999886*aux.LU,0,-0.192755627370783*aux.VU,0]; % Halo
% x0_MCR = [0.191872994233853*aux.LU,0,0,0,-0.531458007318788*aux.VU,0];% 会合周期比 2：1 DRO
% 主星星历积分
[xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);

% 副星星历初值
xx_MCR_temp = Propagate_EphRotFrame(x0_MCR_target,[0,1],[0,1],aux);
x0_MCR_chaser = xx_MCR_temp(2,:);

% x0_TC_VVLH_chaser = T_MCR2VVLH(x0_MCR_chaser-x0_MCR_target,x0_MCR_target);
x0_TC_VVLH_chaser = [1,0.1,0.005,1e-6,1e-6,-1e-06];
% x0_TC_VVLH_chaser = [10,10,10,1e-6,1e-6,-1e-06];
x0_MCR_chaser2 = T_TCO2TCR_eph(x0_TC_VVLH_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;

% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser2,tspan_sec,t_sample,aux);

% 计算相对运动
rvTC_MCR = xx_MCR_chaser-xx_MCR_target;
rvTC_LVLH = T_TCR2TCO_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'LVLH');
var_tested = rvTC_LVLH;
R = diff(var_tested(2:end,1:3))./var_tested(2:end-1,4:6)/diff(t_sample(1:2));
figure(1)
plot(t_sample(3:end)'./86400,R,'.','LineWidth',1.5); 
xlabel('t [day]')
xlim([0,max(t_sample)./86400]); 
ylim([0.95,1.05])
legend('R_x','R_y','R_z')
set(gca,'FontSize',15);


%% CR3BP积分
x0 = [-0.191872994233853,0,0,0,0.531458007318788,0];% 会合周期比 2：1 DRO
T_DRO = 3.394761292341483;%TU
tspan = [0 T_DRO];
options = odeset('RelTol',1e-13,'AbsTol',1e-16);
fun_CR3BP=@(t,x)eom_abs3b(t,x,aux.mu);%相对运动方程
[t_DRO,xx_DRO]=ode45(fun_CR3BP , tspan, x0, options);

%% 画图
f1 = figure(2);
set(f1,'name','星历DRO MCR')

% figure('color',[1 1 1],'name','星历DRO MCR')
p2 = plot3(xx_MCR_target(:,1),xx_MCR_target(:,2),xx_MCR_target(:,3),'color',[0 0.4470 0.7410]); hold on;
plot3(xx_MCR_target(1,1),xx_MCR_target(1,2),xx_MCR_target(1,3),'g^');
plot3(xx_MCR_target(end,1),xx_MCR_target(end,2),xx_MCR_target(end,3),'rv');
% p1 = plot(xx_DRO(:,1)*aux.LU,xx_DRO(:,2)*aux.LU,'k','LineWidth',2); hold on;
% plot(xx_DRO(1,1)*aux.LU,xx_DRO(1,2)*aux.LU,'g^');
box on; grid on; grid minor; hold off;
axis equal; xlabel('x[km]'); ylabel('y[km]')
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
% L = legend([p1 p2],'CR3BP','{星历}DRO');
% set(L,'box','off')
set(gcf,'Renderer','painters');
title('DRO (M frame)')

%% 画图
f2 = figure(3);
set(f2,'name','星历DROrel')
subplot(1,2,1)
plot(rvTC_MCR(:,1),rvTC_MCR(:,2),'color',[0 0.4470 0.7410]); hold on;
plot(rvTC_MCR(1,1),rvTC_MCR(1,2),'g^');
plot(rvTC_MCR(end,1),rvTC_MCR(end,2),'rv');
box on;grid on; grid minor; hold off;
axis equal; xlabel('x[km]'); ylabel('y[km]')
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title('RelMotion(MCR)')

subplot(1,2,2)
plot(rvTC_LVLH(:,1),rvTC_LVLH(:,2),'color',[0 0.4470 0.7410]); hold on;
plot(rvTC_LVLH(1,1),rvTC_LVLH(1,2),'g^');
plot(rvTC_LVLH(end,1),rvTC_LVLH(end,2),'rv');
box on;grid on; grid minor; hold off;
axis equal; xlabel('x[km]'); ylabel('y[km]')
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title('RelMotion(LVLH)')

