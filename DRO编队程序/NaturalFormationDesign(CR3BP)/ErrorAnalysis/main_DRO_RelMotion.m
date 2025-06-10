%% 各个坐标系下的DRO线性化相对运动
% 2022-1-12
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../../subF_eom(CR3BP)')

%% 常数与变量
load('FloquetEig12')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 相对运动积分
x0_REL = Sol_linear.vec2'*1e-3;

% r0_REL = 1e-5*imag(EigenVe(1:3,1)') + 1e-5*k*imag(EigenVe(1:3,5)');
% v0_REL = 1e-5*imag(EigenVe(4:6,1)') + 1e-5*k*imag(EigenVe(4:6,5)');
% x0_REL = [r0_REL,v0_REL];

% dt = 100*para.T0;
dt = 3*para.T0; 
length_t = 20000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;

% 积分
sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
sol_sample = deval(sol,t_sample);
rel_motion_L_linear = sol_sample(7:12,:);
% M坐标系中的DRO绝对运动
abs_motion_M = sol_sample(1:6,:);
% DRO绝对运动转换到惯性系
abs_motion_inertial = synodic2inertial([abs_motion_M(1:3,:)-[1-con.mu;0;0]; abs_motion_M(4:6,:)],t_sample);

% 相对运动坐标转换
rel_motion_M_linear = T_TCR2TCO_CR3BP(rel_motion_L_linear',abs_motion_M','LVLH',con.mu)'; % L→M
rel_motion_VNC_linear = T_TCO2TCR_CR3BP(rel_motion_M_linear',abs_motion_M','VNC',con.mu)'; % M→VNC
rel_motion_VVLH_linear = T_TCO2TCR_CR3BP(rel_motion_M_linear',abs_motion_M','VVLH',con.mu)'; % M→VVLH

%% 画M坐标系和惯性系下的DRO
figure(2)
x_ratio = 3.2; % x轴显示与图形中点的比例
y2x_ratio = 0.5; % 画图时y轴显示与x轴显示的比例
subplot(1,2,1); hold off
plot(abs_motion_inertial(1,1),abs_motion_inertial(2,1),'c*'); hold on
plot(abs_motion_inertial(1,:),abs_motion_inertial(2,:)); hold off
xlabel('x_I'); ylabel('y_I');
legend('InitialPos')
axis equal; title('DRO-Leader (Inertial)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(abs_motion_inertial(1,:));
x_min = min(abs_motion_inertial(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(abs_motion_inertial(2,:))+min(abs_motion_inertial(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

subplot(1,2,2); hold off
plot(abs_motion_M(1,1),abs_motion_M(2,1),'c*'); hold on
plot(abs_motion_M(1,:),abs_motion_M(2,:)); hold off
xlabel('x_M'); ylabel('y_M');
legend('InitialPos')
axis equal; title('DRO-Leader (M)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(abs_motion_M(1,:));
x_min = min(abs_motion_M(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(abs_motion_M(2,:))+min(abs_motion_M(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

%% 各个坐标系下的相对运动
figure(3)
x_ratio = 3; % x轴显示与图形中点的比例
y2x_ratio = 1; % 画图时y轴显示与x轴显示的比例
subplot(2,2,1); hold off
plot(rel_motion_M_linear(1,1),rel_motion_M_linear(2,1),'c*'); hold on
plot(rel_motion_M_linear(1,:),rel_motion_M_linear(2,:)); hold off
xlabel('x_M'); ylabel('y_M');legend('InitialPos')
axis equal; title('RelMotion (M)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(rel_motion_M_linear(1,:));
x_min = min(rel_motion_M_linear(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(rel_motion_M_linear(2,:))+min(rel_motion_M_linear(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

subplot(2,2,2); hold off
plot(rel_motion_L_linear(1,1),rel_motion_L_linear(2,1),'c*'); hold on
plot(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:)); hold off
xlabel('x_L'); ylabel('y_L');legend('InitialPos')
axis equal; title('RelMotion (LVLH)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(rel_motion_L_linear(1,:));
x_min = min(rel_motion_L_linear(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(rel_motion_L_linear(2,:))+min(rel_motion_L_linear(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

subplot(2,2,3); hold off
plot(rel_motion_VNC_linear(1,1),rel_motion_VNC_linear(3,1),'c*'); hold on
plot(rel_motion_VNC_linear(1,:),rel_motion_VNC_linear(3,:)); hold off
xlabel('x_{VNC}'); ylabel('z_{VNC}');legend('InitialPos')
axis equal; title('RelMotion (VNC)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(rel_motion_VNC_linear(1,:));
x_min = min(rel_motion_VNC_linear(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(rel_motion_VNC_linear(3,:))+min(rel_motion_VNC_linear(3,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

subplot(2,2,4); hold off
plot(rel_motion_VVLH_linear(1,1),rel_motion_VVLH_linear(3,1),'c*'); hold on
plot(rel_motion_VVLH_linear(1,:),rel_motion_VVLH_linear(3,:)); hold off
xlabel('x_{VVLH}'); ylabel('z_{VVLH}');legend('InitialPos')
axis equal; title('RelMotion (VVLH)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(rel_motion_VVLH_linear(1,:));
x_min = min(rel_motion_VVLH_linear(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(rel_motion_VVLH_linear(3,:))+min(rel_motion_VVLH_linear(3,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

