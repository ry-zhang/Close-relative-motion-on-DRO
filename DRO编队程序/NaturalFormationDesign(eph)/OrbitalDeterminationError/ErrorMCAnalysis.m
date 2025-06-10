clear
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg
format compact
warning off

UseParallel = 0;

%% 初始误差发散分析
aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;

aux.t0UTC  = [2030 1 1 0 0 0]; % 初始历元
% aux.t0UTC  = [2023 1 1 0 0 0]; % 初始历元

aux = initialize(aux); % 初始化
t0UTC = aux.t0UTC;

% 星历积分初值（旋转系）
x0_MCR_target = [72687.2175909459 0 0 0 -0.540957166578942 0]; % 平面轨道
% 存在地月平面法向振动的轨道
% x0_j2k_target = [380224.0673756608739495, 140818.4956395560002420, 42079.2399940182804130, ...
%     -0.5874909506324951, 0.6787788399867067, 0.3426576490741002];% 星历积分初值（j2000）
% x0_MCR_target = T_ECJ2k2Rot(aux.jd0, x0_j2k_target, [0,0,0],aux.C_Mat, 'MCEMR');

% 主星星历积分（旋转系）
FrameFlag = 'VNC';
tspan_sec = [0,24*86400];% 
size_t = 25;
t_sample = linspace(tspan_sec(1),tspan_sec(2),size_t);
t_sample_jd = aux.jd0 + t_sample/86400;
[xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux,UseParallel);
% 副星星历积分（旋转系）
load('../DROforma2month.mat','MaxDist_all','x0_TC_optimize_all');
scale_0 = 10;
[~,num] = min(MaxDist_all); 
x0_TC_optimize = scale_0*x0_TC_optimize_all(num,:);
x0_TCO_chaser = [x0_TC_optimize(1),0,x0_TC_optimize(2),x0_TC_optimize(3)*1e-6,0,x0_TC_optimize(4)*1e-6];
x0_MCR_chaser = T_TCO2TCR_eph(x0_TCO_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
[xx_MCR_chaser,a_MCR_chaser,STM_j2k] = Propagate_EphRotFrame_STM(x0_MCR_chaser,tspan_sec,t_sample,aux,UseParallel);
xx_TCO_chaser = T_TCR2TCO_eph(xx_MCR_chaser-xx_MCR_target,xx_MCR_target,a_MCR_target,FrameFlag);

xx_j2k_chaser = T_Rot2ECJ2k(t_sample_jd,xx_MCR_chaser,aux.C_Mat,'MCEMR');% km,km/s
x0_j2k_chaser = xx_j2k_chaser(1,:);% km,km/s

%% 误差发散

error_pos = [0.5e-2,3e-2,1e-2];
error_vel = [0.5e-7,1e-7,0.5e-7];
sizeMC = 5000;
xx_TCO_chaser_error_all = zeros(size_t,6,sizeMC);
xx_TCO_chaser_per_all = zeros(size_t,6,sizeMC);
parfor jj_loop = 1:sizeMC
    r_perturb = (2*rand(1,6)-1).*[error_pos,error_vel]; % VVLH中的定轨误差
    
    x0_TCO_chaser_per = x0_TCO_chaser+r_perturb;
    x0_MCR_chaser_per = T_TCO2TCR_eph(x0_TCO_chaser_per,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
    x0_j2k_chaser_per = T_Rot2ECJ2k(aux.jd0,x0_MCR_chaser_per,aux.C_Mat,'MCEMR');% km,km/s
    xx_j2k_chaser_per = zeros(size_t,6);
    for ii_loop = 1:size_t
        phi = reshape(STM_j2k(ii_loop,:),6,6);
        xx_j2k_chaser_per(ii_loop,:) = (x0_j2k_chaser_per-x0_j2k_chaser)*phi'+xx_j2k_chaser(ii_loop,:);
    end
    xx_MCR_chaser_per = T_ECJ2k2Rot(t_sample_jd, xx_j2k_chaser_per, zeros(size_t,3),aux.C_Mat, 'MCEMR');

    xx_TC_MCR_chaser_per = xx_MCR_chaser_per-xx_MCR_target;
    xx_TCO_chaser_per = T_TCR2TCO_eph(xx_TC_MCR_chaser_per,xx_MCR_target,a_MCR_target,FrameFlag);
    xx_TCO_chaser_per_all(:,:,jj_loop) = xx_TCO_chaser_per;
    xx_TCO_chaser_error = xx_TCO_chaser_per - xx_TCO_chaser;
    xx_TCO_chaser_error_all(:,:,jj_loop) = xx_TCO_chaser_error;
end

%% 画误差发散图
% 尝试画包络,调整透明度
figure(10) 
hold off

% for kk_loop = 1:4:size_t
kk_loop_all = 23:1:24;
for kk_loop = kk_loop_all
    xx_TCO_chaser_per_t = reshape(xx_TCO_chaser_per_all(kk_loop,1:3,:),3,sizeMC);
    p1 = scatter3(xx_TCO_chaser_per_t(1,:), xx_TCO_chaser_per_t(2,:), xx_TCO_chaser_per_t(3,:),'.',...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
%     p1.MarkerFaceAlpha = 0.1;
    hold on
end

tspan_sec2 = tspan_sec;% 
size_t2 = 200;
% t_sample2 = linspace(tspan_sec2(1),tspan_sec2(2),size_t2);
t_sample2 = linspace((kk_loop_all(1)-1)*86400,(kk_loop_all(end)-1)*86400,size_t2);
[xx_MCR_target2,a_MCR_target2] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec2,t_sample2,aux,UseParallel);
[xx_MCR_chaser2,a_MCR_chaser2] = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec2,t_sample2,aux,UseParallel);
xx_TCO_chaser2 = T_TCR2TCO_eph(xx_MCR_chaser2-xx_MCR_target2,xx_MCR_target2,a_MCR_target2,FrameFlag);

p2 = plot3(xx_TCO_chaser2(:,1), xx_TCO_chaser2(:,2), xx_TCO_chaser2(:,3),...
    'Color', [200,200,200]/255,'LineWidth',2); hold on;
box on; grid on; grid minor; 

hold off; 
axis equal; xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]')
set(gca,'FontSize',15,'fontname','times new roman'); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title(['RelMotion(',FrameFlag,')'])
% legend1 = legend('0 hrs', '12 hrs','24 hrs','36 hrs','48 hrs');
% legend1 = legend('0 day', '1 day','2 day','3 day','4 day','5 day','6 day');
% legend1 = legend('6 day', '7 day','8 day','9 day','10 day','11 day','12 day');
% legend1 = legend('12 day', '13 day','14 day','15 day','16 day','17 day','18 day');
% legend1 = legend('18 day','19 day', '20 day','21 day','22 day','23 day','24 day');
% legend1.Position = [0.73,0.56,0.16,0.37];
% legendmarkeradjust(12, 'vertical')

% xlim([-22,-12]); zlim([-8,5]);
% view(-40,18)
% view(-8,13)
view([0,0])
% set(gcf,'Renderer','painters')

%% 画时间与速度误差在三个方向的发散极值随时间变化图
figure(11) 
error_max = max(abs(xx_TCO_chaser_error_all),[],3);

subplot(2,1,1)
plot(t_sample/3600, error_max(:,1:3),'LineWidth',2);
xlim([min(t_sample/3600),max(t_sample/3600)])
xlabel('t [hr]'); ylabel('pos error [km]');
set(gca,'FontSize',15,'fontname','times new roman');
title('Evolution of max pos error')
grid on; grid minor
legend('r_x', 'r_y','r_z','Location','northeastoutside');

subplot(2,1,2)
plot(t_sample/3600, error_max(:,4:6),'LineWidth',2);
xlim([min(t_sample/3600),max(t_sample/3600)])
xlabel('t [hr]'); ylabel('vel error [km/s]');
set(gca,'FontSize',15,'fontname','times new roman');
title('Evolution of max vel error')
grid on; grid minor
legend('v_x', 'v_y','v_z','Location','northeastoutside');
