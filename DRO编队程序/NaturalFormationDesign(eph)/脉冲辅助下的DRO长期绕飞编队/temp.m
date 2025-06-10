function [x0f_MCR_target_maneuver,x0f_TCO_rel_maneuver,dv_TCO_all_km,xx_MCR_target_traj,rr_TCO_rel_traj,r_error] ...
    = forcedRelMotion_eph(x0_MCR_target,x_TCO_rel_km,dt_all,aux,Flag_coordinate,isplot,isdisplay)
% 星历模型下DRO附近的受控绕飞
% 
% [x0f_MCR_target_maneuver,x0f_TCO_rel_maneuver,dv_all_km,xx_MCR_target_traj,rr_TCO_rel_traj] = 
% forcedRelMotion_eph(x0_MCR_target,r_TCO_rel_km,dt_all,aux,Flag_coordinate,isplot,isdisplay)
% 
% Input arguments:
% -------------------------------------------------------------------------
% x0_MCR_target         [1x6]       初始时刻的主星状态(km)
% r_TCO_rel_km          [nx3]       副星的变轨点的状态(km) [r0;r1;...;rn](位置/位置+速度均可)
% dt_all                [(n-1)x1]   两次变轨之间的时间间隔(sec)
% Flag_coordinate       [char]      指定副星变轨点的坐标系('LVLH','VNC','VVLH')
% isplot                [1x1]       画图标志位，1是画图，否则不画，默认为0
% isdisplay             [1x1]       显示优化过程标志位，1是显示优化过程，否则不显示，默认为1
% 
% Output arguments:
% -------------------------------------------------------------------------
% x0f_MCR_target_maneuver     [nx6]       变轨点处的主星状态(km)
% x0f_TCO_rel_maneuver        [nx6]       变轨点处变轨后的副星状态(km,km/s)
% dv_TCO_all_km               [nx3]       变轨点处的速度增量(km/s)
% xx_MCR_target_traj                      转移过程中的主星轨迹
% rr_TCO_rel_traj                         转移过程中的相对运动轨迹
% r_error                     [(n-1)x1]   收敛误差
% 
% External functions called:
% -------------------------------------------------------------------------
% maneuver_relmotion, T_TCR2TCO_eph
% 
% Copyright (C) 24/5/2021 by Chihang Yang 
% email: ychhtl@foxmail.com

if nargin == 5
    isplot = 0;
elseif nargin== 6
    isdisplay = 1;
end

%% 常数与变量
opts_fsolve = optimoptions('fsolve','Display','off',...
    'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'Algorithm','levenberg-marquardt','UseParallel',false,...
    'SpecifyObjectiveGradient', true);

%% 正式计算
if length(dt_all) ~= size(x_TCO_rel_km)-1
    error('The number of dt_all should be equal with the number of transfer arcs, as well as the size of r_chaser - 1')
end

% r_TCO_rel_km = r_TCO_rel_km(1:end,:);
t_all = zeros(1,length(dt_all)+1);
for ii_index = 2:length(t_all)
    t_all(ii_index) = sum(dt_all(1:ii_index-1));
end

% 计算变轨处的DRO状态
tspan_sec = [t_all(1),t_all(end)];
[x0f_MCR_target_maneuver,a0f_MCR_target_maneuver] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_all,aux);

%% 将变轨处的目标点转移至VVLH坐标系
if strcmp(Flag_coordinate,'VVLH')
    x_VVLH_rel_all = x_TCO_rel_km;
elseif strcmp(Flag_coordinate,'VNC') || strcmp(Flag_coordinate,'LVLH')
    x_MCR_rel_all = T_TCO2TCR_eph(x_TCO_rel_km,x0f_MCR_target_maneuver,a0f_MCR_target_maneuver,Flag_coordinate);
    x_VVLH_rel_all = T_TCR2TCO_eph(x_MCR_rel_all,x0f_MCR_target_maneuver,a0f_MCR_target_maneuver,'VVLH');
else
    error('Wrong frame');
end

if size(x_VVLH_rel_all,2)==6
    x0_VVLH_rel_loop = x_VVLH_rel_all(1,1:6);
else
    x0_VVLH_rel_loop = [x_VVLH_rel_all(1,1:3),0,0,0];
end


num_dv = size(x_VVLH_rel_all,1)-1;
dv_VVLH_all_km = zeros(num_dv,3);
fval_all = zeros(num_dv,1);
x0f_VVLH_rel_maneuver = zeros(num_dv+1,6);
xx_VVLH_rel_traj = [];
xx_MCR_chaser_traj = [];
xx_MCR_target_traj = [];
aa_MCR_target_traj = [];
if isdisplay == 1
    fprintf('------- 星历下相对运动变轨优化 ---------\n')
end
for ii_index = 1:num_dv
    if isdisplay == 1
        fprintf('第%d/%d段优化 --> ',ii_index,num_dv)
    end
    % 优化
    dt = dt_all(ii_index);
    % 设定初值
    if size(x_VVLH_rel_all,2)==6
        x0 = x_VVLH_rel_all(ii_index,4:6)-x0_VVLH_rel_loop(4:6);
    else
        x0 = [0,0,0];
    end
    rf_VVLH_rel = x_VVLH_rel_all(ii_index+1,1:3);
    x0f_MCR_target = x0f_MCR_target_maneuver([ii_index,ii_index+1],:);
    a0f_MCR_target = a0f_MCR_target_maneuver([ii_index,ii_index+1],:);
    x0_MCR_target_loop = x0f_MCR_target(1,:);
    [dv,fval] = fsolve(@(x)maneuver_relmotion_eph_STM(x,dt,x0f_MCR_target,a0f_MCR_target,x0_VVLH_rel_loop,rf_VVLH_rel,aux),x0,opts_fsolve);

%     paraIter = struct('MaxError',1e-6,'MinDx',1e-15,'BreakError',100,'MaxIter',10,'Display',0);    
%     [dv,fval] = fsolve_GN(@maneuver_relmotion_eph_STM,x0,paraIter,dt,x0f_MCR_target,a0f_MCR_target,x0_REL_loop,rf_REL,aux);

    x0_VVLH_rel_loop(4:6) = x0_VVLH_rel_loop(4:6)+dv;

    % chaser初值, VVLH 2 MCR
    x0_MCR_chaser = T_TCO2TCR_eph(x0_VVLH_rel_loop,x0_MCR_target_loop,a0f_MCR_target(1,:),'VVLH')+x0_MCR_target_loop;
    
    % 副星星历积分
    t_sample = linspace(0,dt,max(ceil(dt/3600)+1,20)); % 一个小时打一个点
    xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,[0 dt],t_sample,aux);
    [xx_MCR_target_loop,aa_MCR_target_loop] = Propagate_EphRotFrame(x0_MCR_target_loop,[0 dt],t_sample,aux);
    % chaser, MCR 2 VVLH
    xx_MCR_rel = xx_MCR_chaser-xx_MCR_target_loop;
    xx_VVLH_rel = T_TCR2TCO_eph(xx_MCR_rel,xx_MCR_target_loop,aa_MCR_target_loop,'VVLH');
    
    if norm(fval)>1e-3
        if nargout < 6
            error('优化未收敛，差%.1e km未到达目标点',norm(fval))
        end
    else
        if isdisplay == 1
            fprintf('优化收敛，收敛误差%.1e km\n',norm(fval))
        end
    end
    
    % 存储数据
    dv_VVLH_all_km(ii_index,:) = dv;
    fval_all(ii_index) = norm(fval);
    x0f_VVLH_rel_maneuver(ii_index,:) = x0_VVLH_rel_loop;
    if ii_index == num_dv
        x0f_VVLH_rel_maneuver(ii_index+1,:) = xx_VVLH_rel(end,1:6);
    end
    
    xx_MCR_chaser_traj = [xx_MCR_chaser_traj;xx_MCR_chaser];
    xx_MCR_target_traj = [xx_MCR_target_traj;xx_MCR_target_loop];
    aa_MCR_target_traj = [aa_MCR_target_traj;aa_MCR_target_loop];
    xx_VVLH_rel_traj = [xx_VVLH_rel_traj;xx_VVLH_rel];
    
    % 更新下一次变轨数据
    x0_VVLH_rel_loop = xx_VVLH_rel(end,1:6);
%     disp([num2str(ii_index),'/',num2str(num_dv)])
    aux.t0UTC = []; % 删除t0UTC以防止与jd0冲突
    aux.jd0 = aux.jd0 + dt / 86400;
end
r_error = fval_all;
% if any(fval_all>1e-3)
%     disp(fval_all);
%     error('The chaser can not arrive the final position(manuver does not converge')
% end
    
if strcmp(Flag_coordinate,'VNC') || strcmp(Flag_coordinate,'LVLH')
    x0f_MCR_rel_maneuver = T_TCO2TCR_eph(x0f_VVLH_rel_maneuver,x0f_MCR_target_maneuver,a0f_MCR_target_maneuver,'VVLH');
    x0f_TCO_rel_maneuver = T_TCR2TCO_eph(x0f_MCR_rel_maneuver,x0f_MCR_target_maneuver,a0f_MCR_target_maneuver,Flag_coordinate);
    xx_MCR_rel_traj = T_TCO2TCR_eph(xx_VVLH_rel_traj,xx_MCR_target_traj,aa_MCR_target_traj,'VVLH');
    rr_TCO_rel_traj = T_TCR2TCO_eph(xx_MCR_rel_traj,xx_MCR_target_traj,aa_MCR_target_traj,Flag_coordinate);
    dv_MCR_all_km = T_TCO2TCR_eph(dv_VVLH_all_km,x0f_MCR_target_maneuver,a0f_MCR_target_maneuver,'VVLH');
    dv_TCO_all_km = T_TCR2TCO_eph(dv_MCR_all_km,x0f_MCR_target_maneuver,a0f_MCR_target_maneuver,Flag_coordinate);
else
    x0f_TCO_rel_maneuver = x0f_VVLH_rel_maneuver;
    rr_TCO_rel_traj = xx_VVLH_rel_traj;
    dv_TCO_all_km = dv_VVLH_all_km;
end


%% 画图
if isplot == 1
    f1 = figure(1);
    set(f1,'name','星历DRO MCR')
    % figure('color',[1 1 1],'name','星历DRO MCR')
    p1 = plot3(xx_MCR_target_traj(:,1),xx_MCR_target_traj(:,2),xx_MCR_target_traj(:,3)); hold on; % ,'color',[0 0.4470 0.7410]
    p2 = plot3(xx_MCR_chaser_traj(:,1),xx_MCR_chaser_traj(:,2),xx_MCR_chaser_traj(:,3));
    p3 = plot3(xx_MCR_target_traj(1,1),xx_MCR_target_traj(1,2),xx_MCR_target_traj(1,3),'g^');
    p4 = plot3(xx_MCR_target_traj(end,1),xx_MCR_target_traj(end,2),xx_MCR_target_traj(end,3),'rv');
    
    box on; grid on; grid minor; hold off;
    axis equal; xlabel('\itx_M \rm[km]'); ylabel('\ity_M \rm[km]'); zlabel('\itz_M \rm[km]')
    L = legend([p1 p2 p3 p4],'target','chaser','Initial Position','Final Position');
%     set(L,'box','off')

    f2 = figure(2);
%     hold off
    p1 = plot3(0,0,0,'ks');  hold on
    p2 = plot3(x0f_TCO_rel_maneuver(1,1),x0f_TCO_rel_maneuver(1,2),x0f_TCO_rel_maneuver(1,3),'g^');
    p3 = plot3(x0f_TCO_rel_maneuver(2:num_dv,1),x0f_TCO_rel_maneuver(2:num_dv,2),x0f_TCO_rel_maneuver(2:num_dv,3),'bo'); 
    p4 = plot3(x0f_TCO_rel_maneuver(end,1),x0f_TCO_rel_maneuver(end,2),x0f_TCO_rel_maneuver(end,3),'rv'); 
    p5 = plot3(rr_TCO_rel_traj(:,1),rr_TCO_rel_traj(:,2),rr_TCO_rel_traj(:,3),'Color',[0, 114, 189]/255,'LineWidth',1.5); 
    hold off
    
%     xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
    
    xlabel('\itx \rm[km]','HorizontalAlignment','center','VerticalAlignment','middle');
    ylabel('\ity \rm[km]','HorizontalAlignment','center','VerticalAlignment','middle'); 
    zlabel('\itz \rm[km]','HorizontalAlignment','center','VerticalAlignment','middle');
    if size(x0f_TCO_rel_maneuver,1) == 2
%         legend('Target S/C','Initial Position','Final Position','Forced Trajectory')
        legend([p1,p2,p4,p5],{'主星','初始点',...
            '目标点','转移轨迹'},'location','eastoutside')
    else
%         legend('Target S/C','Initial Position','Patch Position','Final Position','Forced Trajectory')
        legend([p1,p2,p3,p4,p5],{'主星','初始点',...
            '变轨点','目标点','转移轨迹'},'location','eastoutside')
    end
    axis equal; title(['Forced Motion (',Flag_coordinate,')']); 
    grid on; grid minor
    box on
    % view(0,90)
    % view(38,26)
    f2.Renderer = 'painters';
end
end