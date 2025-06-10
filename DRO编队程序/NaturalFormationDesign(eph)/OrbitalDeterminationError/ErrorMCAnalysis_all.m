clear
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg
format compact
warning off

UseParallel = 1;

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
tspan_sec = [0,6*86400];% 
size_t = 25;
t_sample = linspace(tspan_sec(1),tspan_sec(2),size_t);
t_sample_jd = aux.jd0 + t_sample/86400;
[xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux,UseParallel);
% 副星星历积分（旋转系）
load('../DROforma2month.mat','MaxDist_all','x0_TC_optimize_all');

scale_all = [-100:2:100];
size_scale = length(scale_all);
error_max_all = zeros(size_t,6,size_scale);
for kk_loop = 1:size_scale
    scale_0 = scale_all(kk_loop);
    [~,num] = min(MaxDist_all); 
    x0_TC_optimize = scale_0*x0_TC_optimize_all(num,:);
    x0_TCR_chaser = [x0_TC_optimize(1),0,x0_TC_optimize(2),x0_TC_optimize(3)*1e-6,0,x0_TC_optimize(4)*1e-6];
    x0_MCR_chaser = T_TCO2TCR_eph(x0_TCR_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
    [xx_MCR_chaser,a_MCR_chaser,STM_j2k] = Propagate_EphRotFrame_STM(x0_MCR_chaser,tspan_sec,t_sample,aux,UseParallel);
    xx_TCR_chaser = T_TCR2TCO_eph(xx_MCR_chaser-xx_MCR_target,xx_MCR_target,a_MCR_target,FrameFlag);

    xx_j2k_chaser = T_Rot2ECJ2k(t_sample_jd,xx_MCR_chaser,aux.C_Mat,'MCEMR');% km,km/s
    x0_j2k_chaser = xx_j2k_chaser(1,:);% km,km/s

    %% 误差发散

    error_pos = [0.5e-2,3e-2,1e-2];
    error_vel = [0.5e-7,1e-7,0.5e-7];
    sizeMC = 5000;
    xx_TCR_chaser_error_all = zeros(size_t,6,sizeMC);
    xx_TCR_chaser_per_all = zeros(size_t,6,sizeMC);
    parfor jj_loop = 1:sizeMC
        r_perturb = (2*rand(1,6)-1).*[error_pos,error_vel]; % VVLH中的定轨误差

        x0_TCR_chaser_per = x0_TCR_chaser+r_perturb;
        x0_MCR_chaser_per = T_TCO2TCR_eph(x0_TCR_chaser_per,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
        x0_j2k_chaser_per = T_Rot2ECJ2k(aux.jd0,x0_MCR_chaser_per,aux.C_Mat,'MCEMR');% km,km/s
        xx_j2k_chaser_per = zeros(size_t,6);
        for ii_loop = 1:size_t
            phi = reshape(STM_j2k(ii_loop,:),6,6);
            xx_j2k_chaser_per(ii_loop,:) = (x0_j2k_chaser_per-x0_j2k_chaser)*phi'+xx_j2k_chaser(ii_loop,:);
        end
        xx_MCR_chaser_per = T_ECJ2k2Rot(t_sample_jd, xx_j2k_chaser_per, zeros(size_t,3),aux.C_Mat, 'MCEMR');

        xx_TC_MCR_chaser_per = xx_MCR_chaser_per-xx_MCR_target;
        xx_TCR_chaser_per = T_TCR2TCO_eph(xx_TC_MCR_chaser_per,xx_MCR_target,a_MCR_target,FrameFlag);
        xx_TCR_chaser_per_all(:,:,jj_loop) = xx_TCR_chaser_per;
        xx_TCR_chaser_error = xx_TCR_chaser_per - xx_TCR_chaser;
        xx_TCR_chaser_error_all(:,:,jj_loop) = xx_TCR_chaser_error;
    end
    error_max = max(abs(xx_TCR_chaser_error_all),[],3);
    error_max_all(:,:,kk_loop) = error_max;
    disp([kk_loop size_scale])
end

%% 画48hrs误差随相对运动尺度的关系图
figure(11) 

subplot(2,1,1)
plot(scale_all,reshape(error_max_all(end,1:3,:),3,size_scale),'LineWidth',2);
xlabel('scale'); ylabel('pos error [km]');
set(gca,'FontSize',15,'fontname','times new roman');
ylim([0,0.1])
title('Max pos error at 6 days')
legend('r_x', 'r_y','r_z','Location','northeastoutside');

subplot(2,1,2)
plot(scale_all,reshape(error_max_all(end,4:6,:),3,size_scale),'LineWidth',2);
xlabel('scale'); ylabel('vel error [km/s]');
ylim([0,1.5*max(max(error_max_all(end,4:6,:)))])
set(gca,'FontSize',15,'fontname','times new roman');
title('Max vel error at 6 days')
legend('v_x', 'v_y','v_z','Location','northeastoutside');
