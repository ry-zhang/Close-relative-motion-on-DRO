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
% load('DRO_2to1')
% load('DRO_8toPI')
load('DRO_3to1')
% load('fitted3to1DRO.mat')
% load('DRO_4to1')
x0_DRO_S_2d = state_ini_all(1,:);
x0_DRO_M_3d = [(x0_DRO_S_2d(1:2)-[1-mu,0]),0,x0_DRO_S_2d(3:4),0];
T0 = J_period_all(1,2)*2*pi; % 轨道周期

%%
dt = 1*T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
% t_sample = linspace(0,T0,length_t);
t_sample_day = t_sample*T_norma_day;
T_diff_all = 1./[1e5,10,4,3,2.2,2];
% T_diff_all = 1./[5,4.5,4,2.2,2];
% T_diff_all = [1e-7,1e-1,0.15,0.2,0.247,0.248]; % 1/1e5-1/4.048578,一个自交点
% T_diff_all = [linspace(0.248,0.411,5),0.412]; %1/4.048577~1/2.4317，三个自交点
% T_diff_all = [linspace(0.412,0.494,5),0.495]; % 1/2.4316~1/2.022
% T_diff_all = linspace(0.495,0.5,6); % 1/2.022~1/2
num_Loop = length(T_diff_all);
for kk = 1:num_Loop
    % Leader和follower绝对运动 的积分
    T_diff = T_diff_all(kk);
%     T_diff = 1-0.490;
    sol_abs_l = ode113(@(t,x)eom_abs3b(t,x,mu),[0 dt], x0_DRO_M_3d, opts);% leader abs motion
    xt_DRO_M_3d = deval(sol_abs_l,T0*T_diff);
    sol_abs_f = ode113(@(t,x)eom_abs3b(t,x,mu),[0 dt], xt_DRO_M_3d, opts);% follower abs motion
    sol_abs_l_sample = deval(sol_abs_l,t_sample);
    sol_abs_f_sample = deval(sol_abs_f,t_sample);
    
    % 绝对运动坐标转换
    sol_abs_l_sample_S = sol_abs_l_sample;
    sol_abs_l_sample_S([1,2],:) = sol_abs_l_sample([1,2],:) + [1-mu;0]; % M→S
    abs_motion_inertial = synodic2inertial(sol_abs_l_sample_S,t_sample); % S→I
    % 相对运动坐标转换
    rel_motion_M_nonlinear = sol_abs_f_sample(1:3,:)-sol_abs_l_sample(1:3,:);
    rel_motion_L_nonlinear = T_TCO2TCR_CR3BP(rel_motion_M_nonlinear',sol_abs_l_sample','LVLH',mu)'; % M→L
    [x0,y0,intersectionPos] = selfintersect(rel_motion_L_nonlinear(1,:),rel_motion_L_nonlinear(2,:));
    
   
    %% 转移矩阵的分量、M以及L坐标系中的相对运动分量
    figure(2)
    x_ratio = 3.2; % x轴显示与图形中点的比例
    y2x_ratio = 0.5; % 画图时y轴显示与x轴显示的比例
    subplot(2,2,1); hold off
    plot(abs_motion_inertial(1,1),abs_motion_inertial(2,1),'c*'); hold on
    plot(abs_motion_inertial(1,intersectionPos(:,1)),abs_motion_inertial(2,intersectionPos(:,1)),'bv');
    plot(abs_motion_inertial(1,:),abs_motion_inertial(2,:)); hold off
    xlabel('x_I'); ylabel('y_I');
    legend('InitialPos','IntersecPos')
    axis equal; title('DRO-Leader (Inertial)'); set(gca,'FontSize',13)
    grid on; grid minor
    x_max = max(abs_motion_inertial(1,:));
    x_min = min(abs_motion_inertial(1,:));
    x_middle = (x_max+x_min)/2;
    x_diff = x_max - x_min;
    y_middle = (max(abs_motion_inertial(2,:))+min(abs_motion_inertial(2,:)))/2;
    xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
    ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

    subplot(2,2,2); hold off
    plot(sol_abs_l_sample(1,1),sol_abs_l_sample(2,1),'c*'); hold on
    plot(sol_abs_l_sample(1,intersectionPos(:,1)),sol_abs_l_sample(2,intersectionPos(:,1)),'bv');
    plot(sol_abs_l_sample(1,:),sol_abs_l_sample(2,:)); hold off
    xlabel('x_M'); ylabel('y_M');
    legend('InitialPos','IntersecPos')
    axis equal; title('DRO-Leader (M)'); set(gca,'FontSize',13)
    grid on; grid minor
    x_max = max(sol_abs_l_sample(1,:));
    x_min = min(sol_abs_l_sample(1,:));
    x_middle = (x_max+x_min)/2;
    x_diff = x_max - x_min;
    y_middle = (max(sol_abs_l_sample(2,:))+min(sol_abs_l_sample(2,:)))/2;
    xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
    ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 


    subplot(2,2,3); hold off
    plot(rel_motion_M_nonlinear(1,1),rel_motion_M_nonlinear(2,1),'g^'); hold on
    plot(rel_motion_M_nonlinear(1,intersectionPos(:,1)),rel_motion_M_nonlinear(2,intersectionPos(:,1)),'bv');
    plot(rel_motion_M_nonlinear(1,:),rel_motion_M_nonlinear(2,:)); hold off
    xlabel('x_M'); ylabel('y_M');legend('InitialPos','IntersecPos')
    axis equal; title('RelMotion (M)'); set(gca,'FontSize',13)
    grid on; grid minor
    x_max = max(rel_motion_M_nonlinear(1,:));
    x_min = min(rel_motion_M_nonlinear(1,:));
    x_middle = (x_max+x_min)/2;
    x_diff = x_max - x_min;
    y_middle = (max(rel_motion_M_nonlinear(2,:))+min(rel_motion_M_nonlinear(2,:)))/2;
    xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
    ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

    subplot(2,2,4); hold off
    plot(rel_motion_L_nonlinear(1,1),rel_motion_L_nonlinear(2,1),'ks'); hold on
    plot(rel_motion_L_nonlinear(1,intersectionPos(:,1)),rel_motion_L_nonlinear(2,intersectionPos(:,1)),'bv');
    plot(rel_motion_L_nonlinear(1,:),rel_motion_L_nonlinear(2,:)); hold off
    xlabel('x_L'); ylabel('y_L');legend('InitialPos','IntersecPos')
    axis equal; title('RelMotion (L)'); set(gca,'FontSize',13)
    grid on; grid minor
    x_max = max(rel_motion_L_nonlinear(1,:));
    x_min = min(rel_motion_L_nonlinear(1,:));
    x_middle = (x_max+x_min)/2;
    x_diff = x_max - x_min;
    y_middle = (max(rel_motion_L_nonlinear(2,:))+min(rel_motion_L_nonlinear(2,:)))/2;
    xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
    ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff])    
    %% L坐标系中的相对运动
    x_ratio = 1.6; % x轴显示与图形中点的比例
    y2x_ratio = 4.3; % 画图时y轴显示与x轴显示的比例
    
    figure(5)
    subplot(1,num_Loop,kk)
    plot(rel_motion_L_nonlinear(2,1),rel_motion_L_nonlinear(1,1),'rs'); hold on
    plot(rel_motion_L_nonlinear(2,intersectionPos),rel_motion_L_nonlinear(1,intersectionPos),'bv'); 
    plot(rel_motion_L_nonlinear(2,:),rel_motion_L_nonlinear(1,:)); 
    hold off
    
    if kk==1
        ylabel('x_L');
        legend('InitialPos','IntersecPos') %,'Location','Northeastoutside'
    end
    xlabel('y_L'); 
    axis equal; title([num2str(T_diff,3),'T']); set(gca,'FontSize',13)
    grid on; grid minor
    x_max = max(rel_motion_L_nonlinear(2,:));
    x_min = min(rel_motion_L_nonlinear(2,:));
    x_middle = (x_max+x_min)/2;
    x_diff = x_max - x_min;
    y_middle = (max(rel_motion_L_nonlinear(1,:))+min(rel_motion_L_nonlinear(1,:)))/2;
    xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
    ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff])
    
    pause(0.5)
end

%% 强行滑动x轴分量的相位，得到的现实中不存在的相对运动的图形,
% 以说明即使对于两个相同的坐标轴分量，不同的相位得到的结果也是不同的
% figure(6)
% T_coe_all = [16*5,8,4];
% for jj = 1:3
%     T_coee = T_coe_all(jj);
%     k = length_t/5/T_coee;
%     subplot(1,3,jj)
%     plot(rel_motion_L_nonlinear(1,k+1),rel_motion_L_nonlinear(2,1),'ks'); hold on
%     plot(rel_motion_L_nonlinear(1,k+1:end),rel_motion_L_nonlinear(2,1:length_t-k)); hold off
%     xlabel('x_L'); ylabel('z_L');legend('InitialPos')
%     axis equal; title(['1/',num2str(T_coee),'T ahead']); set(gca,'FontSize',13)
%     grid on; grid minor
% end

%% 保存数据