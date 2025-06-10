%% 三体轨道中的相对运动
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
clear
addpath('../../subF_eom(CR3BP)')

%% 常数与变量
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2
% mu = 0.01211; % 
mu = 0.01215; % 20200531
r_norma = 3.84399*10^5; % km
% T_M = 27.321661; % day
T_norma = sqrt(r_norma^3/(mu_E+mu_M)); % s
T_norma_day = T_norma/3600/24;
v_norma = r_norma/T_norma; % km/s
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 积分与计算
% 绝对运动初值
% 初值由共同质心旋转坐标系S下的动力学（陈冠华程序）导出：原点位于共同质心，x轴从地球指向月球，y轴指向月球绕地球旋转方向
% 绝对动力学位于月球旋转坐标系M：原点位于月球，x轴从地球指向月球，y轴指向月球绕地球旋转的方向
% 惯性坐标系：原点位于共同质心，x、y、z轴与初始时刻的M坐标系相同
% 从S至M的位置转换：rM = rS + [1-mu;0]; 速度转换： vM = vS;
load('DRO_2to1')
x0_DRO_S_2d = state_ini_all(1,:);
x0_DRO_M_3d = [(x0_DRO_S_2d(1:2)-[1-mu,0]),0,x0_DRO_S_2d(3:4),0];
T0 = J_period_all(1,2)*2*pi; % 轨道周期

%%
% 相对运动初值

% r0_REL = [1,0,0] / r_norma/4; % km
% v0_REL = 10*[0,0,-0.284392995600633] * 1e-6/v_norma/4; % mm/s
r0_REL = [-0.683824098909726 ,0,0] * 1e-6; % km
v0_REL = [0,0,0.729646902104155] * 1e-6; % mm/s
    
x0_REL = [r0_REL,v0_REL];

% dt = 5*T0; % 积分时间
dt = 1*T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
% t_sample = linspace(0,T0,length_t);
t_sample_day = t_sample*T_norma_day;
T_diff = 2.08e-4;
% T_diff = 1/8;
delta_tL = 0;
delta_tF = delta_tL + T0/2*T_diff;
% Leader和follower绝对运动 的积分
sol_abs_l = ode113(@(t,x)eom_abs3b(t,x,mu),[0 dt], x0_DRO_M_3d, opts);% leader abs motion
xL0_DRO_M_3d = deval(sol_abs_l,delta_tL);
xF0_DRO_M_3d = deval(sol_abs_l,delta_tF);
sol_abs_l = ode113(@(t,x)eom_abs3b(t,x,mu),[0 dt], xL0_DRO_M_3d, opts);% follower abs motion
sol_abs_f = ode113(@(t,x)eom_abs3b(t,x,mu),[0 dt], xF0_DRO_M_3d, opts);% follower abs motion
sol_abs_l_sample = deval(sol_abs_l,t_sample);
sol_abs_f_sample = deval(sol_abs_f,t_sample);

% 绝对运动坐标转换
sol_abs_l_sample_S = sol_abs_l_sample;
sol_abs_l_sample_S([1,2],:) = sol_abs_l_sample([1,2],:) + [1-mu;0]; % M→S
abs_motion_inertial_l = synodic2inertial(sol_abs_l_sample_S,t_sample+delta_tL); % S→I
sol_abs_f_sample_S = sol_abs_f_sample;
sol_abs_f_sample_S([1,2],:) = sol_abs_f_sample([1,2],:) + [1-mu;0]; % M→S
abs_motion_inertial_f = synodic2inertial(sol_abs_f_sample_S,t_sample+delta_tF); % S→I
% 相对运动坐标转换
rel_motion_M_Abs = sol_abs_f_sample(1:3,:)-sol_abs_l_sample(1:3,:);
rel_motion_L_Abs = T_TCO2TCR_CR3BP(rel_motion_M_Abs',sol_abs_l_sample','LVLH',mu)'; % M→L
rel_motion_L_Abs = rel_motion_L_Abs*r_norma;

%% 画绝对运动轨道、相对运动在L及M坐标系中的轨道，及交叉点对应的位置
% 计算self-intersection points
[x0,y0,intersectionPos] = selfintersect(rel_motion_L_Abs(1,:),rel_motion_L_Abs(2,:));

figure(2)
x_ratio = 3.2; % x轴显示与图形中点的比例
y2x_ratio = 0.5; % 画图时y轴显示与x轴显示的比例
subplot(2,2,1); hold off
plot(abs_motion_inertial_l(1,1),abs_motion_inertial_l(2,1),'c*'); hold on
plot(abs_motion_inertial_f(1,1),abs_motion_inertial_f(2,1),'ro')
plot(abs_motion_inertial_l(1,intersectionPos),abs_motion_inertial_l(2,intersectionPos),'bv');
plot(abs_motion_inertial_l(1,:),abs_motion_inertial_l(2,:)); hold off
xlabel('x_I'); ylabel('y_I');
legend('InitialPos_L','InitialPos_F','IntersecPos')
axis equal; title('DRO-Leader&Follower (Inertial)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(abs_motion_inertial_l(1,:));
x_min = min(abs_motion_inertial_l(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(abs_motion_inertial_l(2,:))+min(abs_motion_inertial_l(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

subplot(2,2,2); hold off
plot(sol_abs_l_sample(1,1),sol_abs_l_sample(2,1),'c*'); hold on
plot(sol_abs_f_sample(1,1),sol_abs_f_sample(2,1),'ro')
plot(sol_abs_f_sample(1,intersectionPos),sol_abs_f_sample(2,intersectionPos),'bv');
plot(sol_abs_l_sample(1,:),sol_abs_l_sample(2,:)); hold off
xlabel('x_M'); ylabel('y_M');
legend('InitialPos_L','InitialPos_F','IntersecPos')
axis equal; title('DRO-Leader&Follower (M)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(sol_abs_f_sample(1,:));
x_min = min(sol_abs_f_sample(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(sol_abs_f_sample(2,:))+min(sol_abs_f_sample(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 


subplot(2,2,3); hold off
plot(rel_motion_M_Abs(1,1),rel_motion_M_Abs(2,1),'g^'); hold on
plot(rel_motion_M_Abs(1,intersectionPos),rel_motion_M_Abs(2,intersectionPos),'bv');
plot(rel_motion_M_Abs(1,:),rel_motion_M_Abs(2,:)); hold off
xlabel('x_M'); ylabel('y_M');legend('InitialPos','IntersecPos')
axis equal; title('RelMotion (M)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(rel_motion_M_Abs(1,:));
x_min = min(rel_motion_M_Abs(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(rel_motion_M_Abs(2,:))+min(rel_motion_M_Abs(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

subplot(2,2,4); hold off
plot(rel_motion_L_Abs(1,1),rel_motion_L_Abs(2,1),'ks'); hold on
plot(rel_motion_L_Abs(1,intersectionPos),rel_motion_L_Abs(2,intersectionPos),'bv');
plot(rel_motion_L_Abs(1,:),rel_motion_L_Abs(2,:)); hold off
xlabel('x_L'); ylabel('y_L');legend('InitialPos','IntersecPos')
axis equal; title('RelMotion (L)'); set(gca,'FontSize',13)
grid on; grid minor
x_max = max(rel_motion_L_Abs(1,:));
x_min = min(rel_motion_L_Abs(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(rel_motion_L_Abs(2,:))+min(rel_motion_L_Abs(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

% subplot(2,2,1); hold on; legend off
% subplot(2,2,2); hold on; legend off
% subplot(2,2,3); hold on; legend off
% subplot(2,2,4); hold on; legend off
% 通过动画表示相对运动在绝对运动轨道中的相对位置

% subplot(2,2,1);h1=line(NaN,NaN,'color','r','linesty','-','LineWidth',1.5);
% h11=line(NaN,NaN,'color','r','Marker','o','linesty','none');
% subplot(2,2,2);h2=line(NaN,NaN,'color','r','linesty','-','LineWidth',1.5);
% h22=line(NaN,NaN,'color','r','Marker','o','linesty','none');
% subplot(2,2,3);h3=line(NaN,NaN,'color','r','linesty','-','LineWidth',1.5);
% h33=line(NaN,NaN,'color','r','Marker','o','linesty','none');
% subplot(2,2,4);h4=line(NaN,NaN,'color','r','linesty','-','LineWidth',1.5);
% h44=line(NaN,NaN,'color','r','Marker','o','linesty','none');

% for anima_index = 1:10:length_t
%     subplot(2,2,1)
%     set(h1,'xdata',abs_motion_inertial_l(1,1:anima_index),'ydata',abs_motion_inertial_l(2,1:anima_index));
%     set(h11,'xdata',abs_motion_inertial_l(1,anima_index),'ydata',abs_motion_inertial_l(2,anima_index));
%     
%     subplot(2,2,2)
%     set(h2,'xdata',sol_abs_l_sample(1,1:anima_index),'ydata',sol_abs_l_sample(2,1:anima_index));
%     set(h22,'xdata',sol_abs_l_sample(1,anima_index),'ydata',sol_abs_l_sample(2,anima_index));
% 
%     subplot(2,2,3)
%     set(h3,'xdata',rel_motion_M_Abs(1,1:anima_index),'ydata',rel_motion_M_Abs(2,1:anima_index));
%     set(h33,'xdata',rel_motion_M_Abs(1,anima_index),'ydata',rel_motion_M_Abs(2,anima_index));
% 
%     subplot(2,2,4)
%     set(h4,'xdata',rel_motion_L_Abs(1,1:anima_index),'ydata',rel_motion_L_Abs(3,1:anima_index));
%     set(h44,'xdata',rel_motion_L_Abs(1,anima_index),'ydata',rel_motion_L_Abs(3,anima_index));
% 
%     pause(0.001)
% end

%% 画图实验，观察L坐标系中，相对运动的x轴分量与z轴分量的关系与联系
figure(11)
subplot(4,1,1)
plot(t_sample/T0,rel_motion_L_Abs(1,:))
ylabel('x_L'); set(gca,'FontSize',13); grid on
subplot(4,1,2)
plot(t_sample/T0,rel_motion_L_Abs(2,:))
ylabel('y_L'); set(gca,'FontSize',13); grid on
subplot(4,1,3)
plot(t_sample/T0,sol_abs_l_sample(1,:))
ylabel('DRO_xM'); set(gca,'FontSize',13); grid on
subplot(4,1,4)
plot(t_sample/T0,sol_abs_l_sample(2,:))
xlabel('t [T]'); ylabel('DRO_yM'); set(gca,'FontSize',13); grid on
