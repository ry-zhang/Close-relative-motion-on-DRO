clear

set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultTextFontSize', 14)
set(0, 'DefaultLineLineWidth', 1)

%% 星历初值设置
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC  = datetime([2030 1 1 0 0 0]); % 初始历元
aux = initialize(aux); % 初始化

% 星历积分区间 
tspan_sec = [0,29.5*86400];% 
t_sample = linspace(tspan_sec(1),tspan_sec(2),500);
% 星历积分初值（旋转系）
x0_MCR_target = [72687.2175909459 0 0 0 -0.540957166578942 0]; % 星历DRO初值
% 主星星历积分
[xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);

% f = dist_DROformation([0.8,0,0,-0.1],xx_MCR_target,a_MCR_target,tspan_sec,t_sample,aux);
% x0_TC_optimize = [0.8,0,0,-0.1];
% lb = []; ub = [];
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point',...
%     'UseParallel',true);
% x0_TC_optimize = fmincon(@(x)dist_DROformation(x,xx_MCR_target,a_MCR_target,tspan_sec,t_sample,aux),...
%     x0,[],[],[],[],lb,ub,[],options);

%% ga优化星历下的自然编队构型
lb = [7,-3,-10,-10]; ub = [13,3,10,10];
options = optimoptions('ga','Display','off','UseParallel',true,'MaxGenerations',100);
x0_TC_optimize_all = zeros(20,4);
MaxDist_all = zeros(20,1);
for ii_loop = 1:20
    tic
    [x0_TC_optimize,fval] = ga(@(x)dist_DROformation(x,xx_MCR_target,a_MCR_target,tspan_sec,t_sample,aux),...
        4,[],[],[],[],lb,ub,[],options);
    x0_TC_optimize_all(ii_loop,:) = x0_TC_optimize;
    MaxDist_all(ii_loop) = fval;
    tf = toc;
    disp([ii_loop,fval,tf])
end
% save DROforma05month MaxDist_all x0_TC_optimize_all



%% 计算旋转系中的运动
load('DROforma2month.mat')
x0_TC_optimize = x0_TC_optimize_all(10,:);
x0_TC_VVLH_chaser = [x0_TC_optimize(1),0,x0_TC_optimize(2),x0_TC_optimize(3)*1e-6,0,x0_TC_optimize(4)*1e-6];
x0_MCR_chaser = T_TCR2TCO_eph(x0_TC_VVLH_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;

% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);

% 计算相对运动
rvTC_MCR = xx_MCR_chaser-xx_MCR_target;
rvTC_VVLH = T_TCO2TCR_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'VVLH');

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
p2 = plot(xx_MCR_target(:,1),xx_MCR_target(:,2),'color',[0 0.4470 0.7410]); hold on;
plot(xx_MCR_target(1,1),xx_MCR_target(1,2),'g^');
plot(xx_MCR_target(end,1),xx_MCR_target(end,2),'rv');
p1 = plot(xx_DRO(:,1)*aux.LU,xx_DRO(:,2)*aux.LU,'k');
box on; grid on; grid minor; hold off;
axis equal; xlabel('x[km]'); ylabel('z[km]')
L = legend([p1 p2],'CR3BP','星历');
set(L,'box','off')

%% 画图
f2 = figure(3);
set(f2,'name','星历DROrel')
subplot(1,2,1)
plot(rvTC_MCR(:,1),rvTC_MCR(:,2),'color',[0 0.4470 0.7410]); hold on;
plot(rvTC_MCR(1,1),rvTC_MCR(1,2),'g^');
plot(rvTC_MCR(end,1),rvTC_MCR(end,2),'rv');
box on;grid on; grid minor; hold off;
axis equal; xlabel('x[km]'); ylabel('y[km]')
title('MCR')
subplot(1,2,2)
plot(rvTC_VVLH(:,1),rvTC_VVLH(:,3),'color',[0 0.4470 0.7410]); hold on;
plot(rvTC_VVLH(1,1),rvTC_VVLH(1,3),'g^');
plot(rvTC_VVLH(end,1),rvTC_VVLH(end,3),'rv');
box on;grid on; grid minor; hold off;
axis equal; xlabel('x[km]'); ylabel('z[km]')
title('VVLH')

