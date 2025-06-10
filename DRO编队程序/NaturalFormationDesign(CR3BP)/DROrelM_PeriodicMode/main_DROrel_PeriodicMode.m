%% 三体轨道中的相对运动
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
clear
addpath('../../subF_eom(CR3BP)')

format longg
format compact
warning off

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

%% 常数与变量
load('FloquetEig12')
% load('FloquetEig12_sy.mat')
% load('generalSolFFT_12')

opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 相对运动积分，计算与坐标转换

% 相对运动初值
% r0_REL = [0.683824098909726 ,0,0] * 1e-6; % 归一化单位
% v0_REL = [0,0,-0.729646902104155] * 1e-6; % 归一化单位

x0_REL = 1e-5*[Sol_linear.vec3]';
% x0_REL = -50/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';
% dt = 5*para.T0; % 积分时间
dt = para.T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;

% 标称轨道及相对运动 的积分
sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
sol_sample = deval(sol,t_sample);
abs_motion_M = sol_sample(1:6,:);
rel_motion_L_linear = sol_sample(7:12,:);

% 绝对运动坐标转换
abs_motion_S = abs_motion_M;
abs_motion_S([1,2],:) = abs_motion_M([1,2],:) + [1-con.mu;0]; % M→S
abs_motion_I = synodic2inertial(abs_motion_S,t_sample); % S→I
% 相对运动坐标转换
rel_motion_M_linear = T_TCR2TCO_CR3BP(rel_motion_L_linear',abs_motion_M','LVLH',con.mu)'; % L→M

%% 画绝对运动轨道、相对运动在L及M坐标系中的轨道，及交叉点对应的位置
%% 计算self-intersection points
% 由于selfintersect函数只能计算小规模的向量，因此叠加做两次交叉点计算
% 第一步，粗算，初步计算出交叉点的位置
[~,~,intersectionPos1] = selfintersect(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:));
% 第二步，细算，在第一步得到的交叉点附近重复迭代计算
t_temple1 = t_sample;
t_temple2 = t_sample;
sol_temple1 = sol_sample;
intersectionPos2 = [intersectionPos1(1),intersectionPos1(2)+length_t];
while norm(sol_temple1(7:9,intersectionPos2(1))-sol_temple1(7:9,intersectionPos2(1)+1))>1e-16
    t_temple1 = linspace(t_temple1(intersectionPos2(1)),t_temple1(intersectionPos2(1)+1),length_t);
    t_temple2 = linspace(t_temple2(intersectionPos2(2)-length_t),t_temple2(intersectionPos2(2)-length_t+1),length_t);
    sol_temple1 = deval(sol,t_temple1);
    sol_temple2 = deval(sol,t_temple2);
    [x0,y0,intersectionPos2] = selfintersect([sol_temple1(7,:),sol_temple2(7,:)],[sol_temple1(8,:),sol_temple2(8,:)]);
end

% 得到两次交叉点对应的时间
t_intersec(1) = t_temple1(intersectionPos2(1));
t_intersec(2) = t_temple2(intersectionPos2(2)-length_t);

DROandRel_intersec = deval(sol,t_intersec);
abs_motion_M_intersec = DROandRel_intersec([1:6],:);
rel_motion_L_linear_intersec = DROandRel_intersec([7:12],:);
% 绝对运动坐标转换
abs_motion_S_intersec = abs_motion_M_intersec;
abs_motion_S_intersec([1,2],:) = abs_motion_M_intersec([1,2],:) + [1-con.mu;0]; % M→S
abs_motion_I_intersec = synodic2inertial(abs_motion_S_intersec,t_intersec); % S→I
% 相对运动坐标转换
rel_motion_M_linear_intersec = T_TCR2TCO_CR3BP(rel_motion_L_linear_intersec',abs_motion_M_intersec','LVLH',con.mu)'; % L→M

%% 将数据单位其转化为km
abs_motion_M_km = abs_motion_M*con.r_norma;
abs_motion_M_intersec_km = abs_motion_M_intersec*con.r_norma;
rel_motion_M_linear_km = rel_motion_M_linear*con.r_norma;
rel_motion_M_linear_intersec_km = rel_motion_M_linear_intersec*con.r_norma;
rel_motion_L_linear_km = rel_motion_L_linear*con.r_norma;
rel_motion_L_linear_intersec_km = rel_motion_L_linear_intersec*con.r_norma;

%% 画图
figure(2)
% x_ratio = 3.2; % x轴显示与图形中点的比例
% y2x_ratio = 0.5; % 画图时y轴显示与x轴显示的比例
x_ratio = 2; y2x_ratio = 1.2;
subplot(1,2,1); hold off
abs_motion_M(1:2,:) = -abs_motion_M(1:2,:);
abs_motion_M_intersec(1:2,:) = -abs_motion_M_intersec(1:2,:);
plot(abs_motion_M(1,1),abs_motion_M(2,1),'g^'); hold on
plot(abs_motion_M_intersec(1,:),abs_motion_M_intersec(2,:),'k*');
plot(abs_motion_M(1,:),abs_motion_M(2,:),'Color',[1, 114, 189]/255,'LineWidth',1.5); hold off
xlabel('\itx_M \rm[LU]'); ylabel('\ity_M \rm[LU]');
legend('InitialPos','IntersecPos')
axis equal; set(gca,'FontSize',13)
grid on; grid minor
x_max = max(abs_motion_M(1,:));
x_min = min(abs_motion_M(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(abs_motion_M(2,:))+min(abs_motion_M(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 
% title('2:1DRO (坐标系M)')
% legend('初始时刻','交叉点')
title('2:1DRO (frame M)')
legend('Initial position','Intersection position')

% subplot(2,2,4); hold off
subplot(1,2,2); hold off
plot(rel_motion_L_linear_km(1,1),rel_motion_L_linear_km(2,1),'g^'); hold on
plot(rel_motion_L_linear_intersec_km(1,:),rel_motion_L_linear_intersec_km(2,:),'k*');
plot(rel_motion_L_linear_km(1,:),rel_motion_L_linear_km(2,:),'Color',[1, 114, 189]/255,'LineWidth',1.5); hold off
xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]');
% legend('InitialPos','IntersecPos')
axis equal; set(gca,'FontSize',13)
grid on; grid minor
x_max = max(rel_motion_L_linear_km(1,:));
x_min = min(rel_motion_L_linear_km(1,:));
x_middle = (x_max+x_min)/2;
x_diff = x_max - x_min;
y_middle = (max(rel_motion_L_linear_km(2,:))+min(rel_motion_L_linear_km(2,:)))/2;
xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 
% title('周期相对运动 (坐标系L)')
% legend('初始时刻','交叉点')
title('periodic mode (frame L)')
legend('Initial position','Intersection position')

% set(gcf,'Color',[255,255,255]/255);
% export_fig PeMode.png -r600

%%
warning on