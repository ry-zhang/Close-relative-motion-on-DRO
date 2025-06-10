%% 三体轨道中的相对运动
% 2021-6-22
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../../subF_eom(CR3BP)')

format longg
format compact

%% 常数与变量
load('FloquetEig12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 相对运动积分，计算与坐标转化
% 相对运动初值
% r0_REL = [0.683824098909726 ,0,0] * 1e-6; % 归一化单位
% v0_REL = [0,0,-0.729646902104155] * 1e-6; % 归一化单位

% x0_REL = -Sol_linear.vec3'*1e-3;
x0_REL = -50/con.r_norma*Sol_linear.vec1'...
    + 50/con.r_norma*Sol_linear.vec5';
% dt = 5*para.T0; % 积分时间
dt = 4*para.T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;

% 标称轨道及相对运动 的积分
sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
sol_sample = deval(sol,t_sample);
abs_motion_M = sol_sample(1:6,:);
rel_motion_L_linear = sol_sample(7:12,:);

rel_motion_L_linear_km = rel_motion_L_linear(1:3,:)*con.r_norma;

%% 画图
figure(2)
plot3(rel_motion_L_linear_km(1,:),rel_motion_L_linear_km(2,:),rel_motion_L_linear_km(3,:),'Color',[237, 177, 32]/255,'LineWidth',1.5); hold on
p1 = plot3(rel_motion_L_linear_km(1,1),rel_motion_L_linear_km(2,1),rel_motion_L_linear_km(3,1),'g^');
p2 = plot3(rel_motion_L_linear_km(1,end),rel_motion_L_linear_km(2,end),rel_motion_L_linear_km(3,end),'rv');
hold off
xlabel('x_L [km]'); ylabel('y_L [km]'); zlabel('z_L [km]');
legend([p1,p2],'InitialPos','FinalPos')
axis equal; title('RelMotion (L frame)'); set(gca,'FontSize',13)
grid on; grid minor

% 各分量图
figure(3)
plot(t_sample_day,rel_motion_L_linear_km,'LineWidth',1.5)
legend('x','y','z')
xlabel('t [day]'); ylabel('components [km]');
set(gca,'FontSize',13)
grid on; grid minor