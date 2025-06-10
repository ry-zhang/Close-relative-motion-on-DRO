function [x0_DRO_all,x0_REL_all_km,dv_all_km,r_REL_traj_all_km,x_DRO_traj_all] = forcedRelMotion(x0_target,r_chaser_km,dt_all,Flag_coordinate,isplot,isDisplay)
% CR3BP模型下DRO附近的受控绕飞
% 
% [x0_DRO_all,x0_REL_all_km,dv_all_km] = 
% forcedRelMotion(x0_target,r_chaser_km,dt_all,Flag_coordinate,isplot)
%
% Input arguments:
% -------------------------------------------------------------------------
% x0_target             [1x6]       初始时刻的主星状态(归一化单位)
% r_chaser_km    [nx3]       副星的变轨点的位置(km) [r0;r1;...;rn]
% dt_all                [(n-1)x1]   两次变轨之间的时间间隔
% Flag_coordinate       [char]      指定副星变轨点的坐标系('LVLH','VNC','VVLH')
% isplot                [1x1]       画图标志位，1是画图，否则不画
% UseParallel           [1x1]       并行标志位，1是并行，否则不并行
% 
% Output arguments:
% -------------------------------------------------------------------------
% x0_DRO_all            [nx6]       变轨点处的主星状态(归一化单位)
% x0_REL_all_km         [nx6]       变轨点处变轨后的副星状态(km,km/s)
% dv_all_km             [nx3]       变轨点处的速度增量(km/s)
% 
% External functions called:
% -------------------------------------------------------------------------
% eom_abs3b, eom_rel3b, maneuver_relmotion, T_TCR2TCO_CR3BP, T_TCO2TCR_CR3BP
% 
% Copyright (C) 24/5/2021 by Chihang Yang 
% email: ychhtl@foxmail.com
% 

if nargin == 4
    isplot = 0;
elseif nargin == 5
    isDisplay = 1;
end
%% 常数与变量
para.mu = 0.01215; % 20200531
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2
para.r_norma = 3.84399*10^5; % km
para.T_norma = sqrt(para.r_norma^3/(mu_E+mu_M)); % s
para.T_norma_day = para.T_norma/3600/24;
para.v_norma = para.r_norma/para.T_norma; % km/s

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);
opts_fsolve = optimoptions('fsolve','Display','off',...
    'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'Algorithm','levenberg-marquardt','UseParallel',false);

%% 正式计算
if length(dt_all) ~= size(r_chaser_km)-1
    error('The number of dt_all should be equal with the number of transfer arcs, as well as the size of r_chaser - 1')
end

r_chaser_all = r_chaser_km(1:end,:)/para.r_norma;
t_all = zeros(1,length(dt_all)+1);
for ii_index = 2:length(t_all)
    t_all(ii_index) = sum(dt_all(1:ii_index-1));
end

% 计算变轨处的DRO状态
[~,x0_DRO_all] = ode113(@(t,x)eom_abs3b(t,x,para.mu),t_all, x0_target, opts_ode);
if size(t_all,2) == 2
    x0_DRO_all = x0_DRO_all([1,end],:);
end

% 将变轨处的目标点转移至LVLH坐标系
if strcmp(Flag_coordinate,'LVLH')  % VNC 2 LVLH
    r_chaser_L_all = r_chaser_all;
elseif strcmp(Flag_coordinate,'VNC') || strcmp(Flag_coordinate,'VVLH')
    r_chaser_M = T_TCO2TCR_CR3BP(r_chaser_all,x0_DRO_all,Flag_coordinate,para.mu);
    r_chaser_L_all = T_TCR2TCO_CR3BP(r_chaser_M,x0_DRO_all,'LVLH',para.mu);
else
    error('Wrong orbital frame flag')
end

x0_DRO_loop = x0_target;
if size(r_chaser_L_all,2) == 3
    x0_REL_loop = [r_chaser_L_all(1,:),0,0,0];
else
    x0_REL_loop = r_chaser_L_all(1,:);
end
num_dv = size(r_chaser_L_all,1)-1;
dv_all_L = zeros(num_dv,3);
fval_all = zeros(num_dv,1);
x0_REL_L_all = zeros(num_dv+1,6);
r_REL_traj_L_all = [];
if isDisplay == 1
    fprintf('------- CRTBP下相对运动变轨优化 ---------\n')
end

for ii_index = 1:num_dv
    if isDisplay == 1
        fprintf('第%d段优化,共%d段优化 --> ',ii_index,num_dv)
    end
    % 优化
    dt = dt_all(ii_index);
    x0 = [0,0,0];
    rf_rel = r_chaser_L_all(ii_index+1,:);
    [dv,fval] = fsolve(@(x)maneuver_relmotion(x,dt,x0_DRO_loop,x0_REL_loop,rf_rel,para),x0,opts_fsolve);
    x0_REL_loop(4:6) = x0_REL_loop(4:6)+dv;

    % 相对运动在LVLH的积分
    sol = ode113(@(t,x)eom_rel3b(t,x,para.mu),[0 dt], [x0_DRO_loop, x0_REL_loop], opts_ode);
    t_sample = linspace(0,dt,max(ceil(dt/(3600/86400/27.28)),1000)); % 一个小时打一个点
    y_temp = deval(sol,t_sample)';
    
    % 存储数据
    dv_all_L(ii_index,:) = dv;
    fval_all(ii_index) = norm(fval);
    x0_REL_L_all(ii_index,:) = x0_REL_loop;
    if ii_index == num_dv
        x0_REL_L_all(ii_index+1,:) = [y_temp(end,7:9),y_temp(end,10:12)];
    end
    r_REL_traj_L_all = [r_REL_traj_L_all;y_temp];
    
    % 更新下一回合数据
    x0_DRO_loop = y_temp(end,1:6);
    x0_REL_loop = y_temp(end,7:12);
    if norm(fval)>1e-5
        error('优化未收敛，差%.1e km未到达目标点',norm(fval))
    else
        if isDisplay == 1
            fprintf('优化收敛，收敛误差%.1e km\n',norm(fval))
        end
    end
end

% if any(fval_all>1e-9)
%     error('The chaser can not arrive the final position(manuver does not converge')
% else
%     disp('------- 优化收敛 ---------')
% end

if strcmp(Flag_coordinate,'LVLH') % VNC 2 LVLH
    r_REL_traj_all = r_REL_traj_L_all(:,7:9);
    x0_REL_all = x0_REL_L_all;
    dv_all = dv_all_L;
elseif strcmp(Flag_coordinate,'VNC') || strcmp(Flag_coordinate,'VVLH')
    x0_REL_M_all = T_TCO2TCR_CR3BP(x0_REL_L_all,x0_DRO_all,'LVLH',para.mu);
    x0_REL_all = T_TCR2TCO_CR3BP(x0_REL_M_all,x0_DRO_all,Flag_coordinate,para.mu);
    r_REL_traj_M_all = T_TCO2TCR_CR3BP(r_REL_traj_L_all(:,7:9),r_REL_traj_L_all(:,1:6),'LVLH',para.mu);
    r_REL_traj_all = T_TCR2TCO_CR3BP(r_REL_traj_M_all,r_REL_traj_L_all(:,1:6),Flag_coordinate,para.mu);
    dv_all_M = T_TCO2TCR_CR3BP(dv_all_L,x0_DRO_all,'LVLH',para.mu);
    dv_all = T_TCR2TCO_CR3BP(dv_all_M,x0_DRO_all,Flag_coordinate,para.mu);
end

dv_all_km = dv_all*para.v_norma
r0_REL_all_km = x0_REL_all(:,1:3)*para.r_norma;
x0_REL_all_km = [r0_REL_all_km,x0_REL_all(:,4:6)*para.v_norma];
r_REL_traj_all_km = r_REL_traj_all*para.r_norma;
x_DRO_traj_all = r_REL_traj_L_all(:,1:6);
% dvnorm_all_km = sqrt(sum(dv_all_km.^2,2));
% dv_total_km = sum(dvnorm_all_km);

% 画图也要判断
%% 画图
if isplot == 1
    fig = figure(1);
%     hold off
    p1 = plot3(0,0,0,'ks');  hold on
    p2 = plot3(r0_REL_all_km(1,1),r0_REL_all_km(1,2),r0_REL_all_km(1,3),'g^');
    p3 = plot3(r0_REL_all_km(2:num_dv,1),r0_REL_all_km(2:num_dv,2),r0_REL_all_km(2:num_dv,3),'bo'); 
    p4 = plot3(r0_REL_all_km(end,1),r0_REL_all_km(end,2),r0_REL_all_km(end,3),'rv'); 
    colorTable = get(gca,'colororder');
    p5 = plot3(r_REL_traj_all_km(:,1),r_REL_traj_all_km(:,2),r_REL_traj_all_km(:,3),'Color',colorTable(1,:),'LineWidth',1.5); 
    hold off
    xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]'); zlabel('\itz \rm[km]');
%     xlabel('\itx \rm[km]','HorizontalAlignment','center','VerticalAlignment','middle');
%     ylabel('\ity \rm[km]','HorizontalAlignment','center','VerticalAlignment','middle'); 
%     zlabel('\itz \rm[km]','HorizontalAlignment','center','VerticalAlignment','middle');
    if size(r0_REL_all_km,1) == 2
%         legend('Target S/C','Initial Position','Final Position','Forced Trajectory')
        legend([p1,p2,p4,p5],{'主星','初始点',...
            '目标点','转移轨迹'},'location','eastoutside')
    else
%         legend('Target S/C','Initial Position','Patch Position','Final Position','Forced Trajectory')
        legend([p1,p2,p3,p4,p5],{'主星','初始点',...
            '变轨点','目标点','转移轨迹'},'location','eastoutside')
    end
    axis equal; title(['Forced Motion (',Flag_coordinate,',CR3BP)']); 
%     set(gca,'FontSize',15)
    grid on; grid minor
    box on

    fig.Renderer = 'painters';
end
end