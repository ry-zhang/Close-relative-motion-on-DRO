%% CR3BP下的强迫相对运动
% 2021-7-11
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../subF_eom(CR3BP)')

format longg
format compact
% warning off
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'DefaultAxesFontsize',15);
set(0,'DefaultTextFontsize',15);
set(0,'defaultLineLineWidth',1.5)

%% 常数与变量

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);

% 加载周期轨道及相对运动轨道
load('FloquetEig12.mat')
x0_DRO_0 = x0_DRO_M_3d;
x0_REL_period = -2/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';
x0_target = x0_DRO_0;

Flag_coordinate = 'LVLH';
% Flag_coordinate = 'VVLH';
% Flag_coordinate = 'VNC';

% isplot = 1;
isplot = 0;

%% 选择patch pointd 的方式
patchflag = 1;
if patchflag == 1 %% circumnavigation，绕飞
    n_samp = 10; % 每圈的点数
    dt_day = 25;
    dt = dt_day/con.T_norma_day/(n_samp-1); % 转移时间
    % 空间圆方程
    radi = 10; % 半径，km
    r_cent = [0,0,0]; % 圆心，km
%     r_normal = [0,0,1]; % 法向量
    r_normal = [0,1,1]; % 法向量
    r_normal = r_normal/norm(r_normal); % 单位法向量
    % 计算r_a与r_b,先对r_normal升序排列，以寻找最大的非零值，预防奇异
    [r_normal_resort,r_normal_index] = sort(abs(r_normal),'ascend');
    r_a1 = 0.1;
    r_a2 = 0.2;
    r_a3 = -(r_a1*r_normal_resort(1)+r_a2*r_normal_resort(2))/r_normal_resort(3);
    r_a_temp = [r_a1,r_a2,r_a3];
    eye3 = eye(3);
    r_a = r_a_temp*eye3(r_normal_index,:)/norm(r_a_temp);
    r_b = cross(r_normal,r_a);
    theta_all = linspace(0,2*pi*n_samp/(n_samp+1),n_samp);
    r_chaser_km = r_cent + radi*cos(theta_all)'*r_a + radi*sin(theta_all)'*r_b;
elseif patchflag == 2   %% 点到点自主分段转移
    r0 = [-3,1,2];
    v0 = [0,0,0];
    rf = [1,2,-2];
    
    dt_day = 6;
    n_samp = 6; % 变轨点数，若点数为2，则是两点间直接转移
    dt = dt_day/con.T_norma_day/(n_samp-1);

    % 计算旋转轴和夹角
    if norm(cross(r0,rf))>=1e-8 % r0与rf方向不同
        r_rot = cross(r0,rf)/norm(cross(r0,rf));
        theta = acos(dot(r0,rf)/norm(r0)/norm(rf));
    else
        [r0_com,order] = sort(abs(r0),'descend');
        if norm(r0)<=1e-8 || norm(rf)<=1e-8 || r0(order(1))/rf(order(1))>0 % r0与rf中存在零向量，或二者同向
            r_rot = [0,0,0];
            theta = pi/4;
        else % r0(order(1))/rf(order(1))<0 % r0与rf反向
            r_rot(order(2:3)) = [1,1];
            r_rot(order(1)) = dot(r_rot(order(2:3)),r0(order(2:3)))/r0(order(1));
            r_rot = r_rot/norm(r_rot);
            theta = acos(dot(r0,rf)/norm(r0)/norm(rf));
        end
    end
    % 将r0绕r_rot旋转至rf，同时进行线性放缩，以得到一条从r0至rf的渐近曲线
    theta_all = linspace(0,theta,n_samp);
    r_chaser_km = zeros(n_samp,3); % 计算拼接点
    for jj_loop = 1:n_samp
        theta_temp = theta_all(jj_loop);
        if norm(r0) > 1e-8 % 将r0 绕r_rot旋转theta_temp
            r_patch_temp = r0*cos(theta_temp)+cross(r_rot,r0)*sin(theta_temp)...
                +dot(r_rot,r0)*r_rot*(1-cos(theta_temp));
        else % 此时r_rot = [0,0,0],从r0直线趋近于rf
            theta_temp2 = theta_temp-theta;
            r_patch_temp = rf*cos(theta_temp2);
        end
        if norm(r_patch_temp) > 1e-8  % 归一化
            r_patch_temp = r_patch_temp./norm(r_patch_temp);
        end
        scale = norm(r0)*(n_samp - jj_loop)/(n_samp-1)+norm(rf)*(jj_loop-1)/(n_samp-1);
        r_chaser_km(jj_loop,:) = r_patch_temp*scale;
    end

else  %% 定点悬停
%     dt = 24/24/con.T_norma_day; % 转移时间
%     r_chaser_km = [3,3,3; ];
end

%% 计算变轨
dt_all = dt*ones(size(r_chaser_km,1)-1,1); % 时间

[x0_DRO_all,r0_REL_all_km,dv_all_km,r_REL_traj_all_km] = forcedRelMotion(x0_target,r_chaser_km,dt_all,Flag_coordinate,isplot,0);
dvnorm_all_km = sqrt(sum(dv_all_km.^2,2));
dv_total_m = sum(dvnorm_all_km)*1e3;

% 画圆型参考轨迹
figure(1)
% hold on
p1 = plot3(0,0,0,'ks');  hold on
p2 = plot3(r0_REL_all_km(1,1),r0_REL_all_km(1,2),r0_REL_all_km(1,3),'g^');
p3 = plot3(r0_REL_all_km(2:end-1,1),r0_REL_all_km(2:end-1,2),r0_REL_all_km(2:end-1,3),'bo'); 
p4 = plot3(r0_REL_all_km(end,1),r0_REL_all_km(end,2),r0_REL_all_km(end,3),'rv'); 
colorTable = get(gca,'colororder');
p5 = plot3(r_REL_traj_all_km(:,1),r_REL_traj_all_km(:,2),r_REL_traj_all_km(:,3),'Color',colorTable(1,:),'LineWidth',1.5); 
% n_samp_ref = 180;
% theta_all = linspace(0,2*pi*n_samp_ref/(n_samp_ref+1),n_samp_ref);
% r_chaser_km_ref = r_cent + radi*cos(theta_all)'*r_a + radi*sin(theta_all)'*r_b;
% p6 = plot3(r_chaser_km_ref(:,1),r_chaser_km_ref(:,2),r_chaser_km_ref(:,3),'LineWidth',1,'Color',[100, 100, 100]/255);

%     plot3(r_chaser_km(:,1),r_chaser_km(:,2),r_chaser_km(:,3),'LineWidth',1.5,'Color',[0, 114, 189]/255); axis equal
plot3([0,r_chaser_km(1,1)],[0,r_chaser_km(1,2)],[0,r_chaser_km(1,3)],'LineWidth',1.5,'Color',[236, 119, 40]/255);
plot3([0,r_chaser_km(end,1)],[0,r_chaser_km(end,2)],[0,r_chaser_km(end,3)],'LineWidth',1.5,'Color',[236, 119, 40]/255);

hold off
view([117,25])
% view([-123,21])
% view([0,90])
label('\itx_L \rm[km]','\ity_L \rm[km]','\itz_L \rm[km]'); title('')
% legend([p1,p2,p3,p4,p5,p6],{'主星','初始点','变轨点','目标点','转移轨迹','参考轨迹'},'location','eastoutside')
legend([p1,p2,p3,p4,p5],{'主星','初始点','变轨点','目标点','转移轨迹'},'location','eastoutside')
axis equal;
grid on; grid minor; box on

% title({'Forced Motion',['(',Flag_coordinate,', CR3BP, ',num2str(dt_day), ' days)']});
title({'受控绕飞编队',['(CRTBP, ', num2str(dt_day), '天)']}); 
% title({'两点安全转移轨迹',['(CRTBP, ', num2str(dt_day), '天)']}); 

