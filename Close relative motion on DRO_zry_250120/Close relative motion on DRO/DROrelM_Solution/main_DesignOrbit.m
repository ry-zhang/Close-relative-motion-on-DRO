close all;
%%设计DRO轨道，6种模态组合
%2024-11-13
% by ZhangRuYue
% email:zhangruyue22@csu.ac.cn
addpath('../subF_eom(CR3BP)')
%% 常数与变量
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2  //4902.56146783419
% con.mu = 0.01211; % 
con.mu = 0.01215; % 20200531
con.r_norma = 3.84399*10^5; % km
%T_M = 27.321661; % day
con.T_norma =2*pi*sqrt(con.r_norma^3/(mu_E+mu_M)); % s
con.T_norma_day = con.T_norma/3600/24;
con.v_norma = con.r_norma/con.T_norma; % km/s
%con.v_norma = con.r_norma/(T_M*24*2600); % km/s
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);
load('DRO_all_MoonPeriod')

numDRO = 4;%2:1
x0_DRO_S_2d = state_ini_all(numDRO,:);
x0_DRO_M_3d = [(x0_DRO_S_2d(1:2)-[1-con.mu,0]),0,x0_DRO_S_2d(3:4),0];
para.T0 = J_period_all(numDRO,2)*2*pi; % 轨道周期

%% 6种单一模态初值
% load(Sol_linear,'Sol_linear.vec');
% load(EigenVector,'EigenVector.vec');
x0_REL= [];
for jj_index = 1:6
        eval(['x0 = Sol_linear.vec',num2str(jj_index),'*1e-5;'])
        x0_REL(jj_index,:) =x0';
end

C1=0;
C2=0;
C3=0;
C4=1;
C5=0;
C6=0;
%%混合模态初值
x0_REL_MIX=[];
x0_REL_MIX(1,:)=C1*x0_REL(1,:)+C2*x0_REL(2,:) +C3*x0_REL(3,:) +C4*x0_REL(4,:) +C5*x0_REL(5,:)+C6*x0_REL(6,:);
N=50/x0_REL_MIX(2)/con.r_norma;
x0_REL_MIX=N*x0_REL_MIX;

%% 相对运动积分

dt = 5*para.T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;

% 周期轨道
% 标称轨道及相对运动 p3 的积分
sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [para.x0_DRO, x0_REL_MIX], opts);
sol_sample = deval(sol,t_sample);
abs_motion_L_M = sol_sample(1:6,:);
rel_motion_L_linear = [sol_sample(7:9,:)*con.r_norma;sol_sample(10:12,:)*con.v_norma];
 % 将绝对运动轨道转移至惯性系中
 abs_motion_M = sol_sample(1:6,:);

%% 画图
isplot=1;
if isplot == 1
    dt = 5*para.T0;
    % dt = 1*T0;
    length_t = 2000;
    t_sample = linspace(0,dt,length_t);
    opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
    % 假设 rel_motion_L_linear 是一个矩阵
    n_cols = size(rel_motion_L_linear, 2);  % 获取 rel_motion_L_linear 第3行的列数

    % 创建一个与 rel_motion_L_linear 第3行列数相同的全1向量
    rel_motion_L_linear_z = rel_motion_L_linear(3,end)*ones(1, n_cols);  % 初始化为 1 的向量，行数为1，列数与 rel_motion_L_linear(3,:) 相同

    figure(1)
    plot3(rel_motion_L_linear(1,1),rel_motion_L_linear(2,1),rel_motion_L_linear(3,1),'g^'); hold on
    plot3(rel_motion_L_linear(1,end),rel_motion_L_linear(2,end),rel_motion_L_linear(3,end),'r*');  hold on
    plot3(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:),rel_motion_L_linear(3,:),'Color',[0, 114, 189]/255,'LineWidth',1.5) ;hold on
    plot3(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:),rel_motion_L_linear_z(1,:),'Color',[0, 114/255, 189/255,0.5],'LineWidth',0.5); hold on

    legend('初始时刻','终端时刻','三维相对轨迹','平面相对轨迹','Location', 'northeast','fontsize',8,'interpreter','latex')
%     legend('初始时刻','终端时刻')
    xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]');zlabel('\itz \rm[km]');
    axis equal;
    x_ratio = 2; % x轴显示与图形中点的比例
    y2x_ratio = 0.7; % 画图时y轴显示与x轴显示的比例
    x_max = max(rel_motion_L_linear(1,:));
    x_min = min(rel_motion_L_linear(1,:));
    x_middle = (x_max+x_min)/2;
    x_diff = x_max - x_min;
    y_middle = (max(rel_motion_L_linear(2,:))+min(rel_motion_L_linear(2,:)))/2;
    xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
    ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 
    title({'月心LVLH下相对轨迹'},'fontsize',16,'interpreter','latex');hold on;
    text(0.9, 0.5, '1*x3+1*x5混合模态', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'interpreter', 'latex');
    

    grid on; grid minor

end

