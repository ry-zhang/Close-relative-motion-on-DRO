%% 三体轨道中的相对运动
% 2021-03-24
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
% load('FloquetEig12')
load('FloquetEig12_sy.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 相对运动积分，计算与坐标转化
% 相对运动初值
% r0_REL = [0.683824098909726 ,0,0] * 1e-6; % 归一化单位
% v0_REL = [0,0,-0.729646902104155] * 1e-6; % 归一化单位
num = 300;
r_ratio_all = linspace(-300,300,num)/con.r_norma/Sol_linear.vec3(2);
T_in_all = zeros(num,1);
T_out_all = zeros(num,1);
x_in_amp_all = zeros(num,1);
y_in_amp_all = zeros(num,1);
x_out_amp_all = zeros(num,1);
y_out_amp_all = zeros(num,1);
dist_min_all = zeros(num,1);
dist_max_all = zeros(num,1);

parfor ii_index = 1:num
    x0_REL = -Sol_linear.vec3'*r_ratio_all(ii_index);
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

    %% 画绝对运动轨道、相对运动在L及M坐标系中的轨道，及交叉点对应的位置
    %% 计算self-intersection points
    % 由于selfintersect函数只能计算小规模的向量，因此叠加做两次交叉点计算
    % 第一步，粗算，初步计算出交叉点的位置
    [~,~,intersectionPos1] = selfintersect(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:));
    t_intersec = zeros(1,2);
    t_intersec(1) = t_sample(intersectionPos1(1));
    t_intersec(2) = t_sample(intersectionPos1(2));
    
    %% 将数据单位其转化为km
    rel_motion_L_linear_km = rel_motion_L_linear*con.r_norma;
    
    %% 相位差，最大最远距离存储
    rel_motion_L_linear_in = rel_motion_L_linear_km(:,intersectionPos1(1):intersectionPos1(2));
    rel_motion_L_linear_out = rel_motion_L_linear_km(:,[1:intersectionPos1(1),intersectionPos1(2):end]);
    dist_all = sqrt(sum(rel_motion_L_linear_km.^2));
    
    T_in_all(ii_index) = diff(t_intersec)*con.T_norma_day;
    T_out_all(ii_index) = (para.T0-diff(t_intersec))*con.T_norma_day;
    x_in_amp_all(ii_index) = max(rel_motion_L_linear_in(1,:))-min(rel_motion_L_linear_in(1,:));
    y_in_amp_all(ii_index) = max(rel_motion_L_linear_in(2,:))-min(rel_motion_L_linear_in(2,:));
    x_out_amp_all(ii_index) = max(rel_motion_L_linear_out(1,:))-min(rel_motion_L_linear_out(1,:));
    y_out_amp_all(ii_index) = max(rel_motion_L_linear_out(2,:))-min(rel_motion_L_linear_out(2,:));
    dist_min_all(ii_index) = min(dist_all);
    dist_max_all(ii_index) = max(dist_all);
end

%% 画图
% 内外圈所需时间
y0_all = r_ratio_all*con.r_norma*Sol_linear.vec3(2);
% y0_all = r_ratio_all;
figure(1); hold off
plot(y0_all,T_in_all,'LineWidth',1.5); hold on
plot(y0_all,T_out_all,'LineWidth',1.5);
xlabel('\ity\rm(0) [km]'); 
% xlabel('\itk_{\rm0}'); 
ylabel('\itt \rm[day]');
% legend('Inner circle','Outer circle')
legend('内圈','外圈')
xlim([min(y0_all),max(y0_all)]); ylim([6.5,8.5])
set(gca,'FontSize',15)
grid on; grid minor

figure(2); hold off
plot(y0_all,dist_min_all,'LineWidth',1.5); hold on
plot(y0_all,dist_max_all,'LineWidth',1.5);
% xlabel('\itk_{\rm0}'); 
xlabel('\ity\rm(0) [km]'); 
ylabel('distance [km]');
% legend('minimal distance','maximal distance','Location','north')
legend('最小距离','最大距离','Location','north')
xlim([min(y0_all),max(y0_all)]); 
set(gca,'FontSize',15)
grid on; grid minor

figure(3); hold off
plot(y0_all,x_in_amp_all,'-.','Color',[1, 114, 189]/255,'LineWidth',1.5); hold on
plot(y0_all,x_out_amp_all,'-.','Color',[217, 83, 25]/255,'LineWidth',1.5);
plot(y0_all,y_in_amp_all,'Color',[1, 114, 189]/255,'LineWidth',1.5);
plot(y0_all,y_out_amp_all,'Color',[217, 83, 25]/255,'LineWidth',1.5);
xlabel('\ity\rm(0) [km]'); 
% xlabel('\itk_{\rm0}'); 
ylabel('距离 [km]');
% legend('x-amplitude (inner circle)','x-amplitude (outer circle)',...
%     'y-amplitude (inner circle)','y-amplitude (outer circle)',...
%     'Location','north')
legend('\itx\rm振幅(内圈)','\itx\rm振幅(外圈)',...
    '\ity\rm振幅(内圈)','\ity\rm振幅(外圈)',...
    'Location','north')
xlim([min(y0_all),max(y0_all)]); 
set(gca,'FontSize',15)
grid on; grid minor

%%
warning on