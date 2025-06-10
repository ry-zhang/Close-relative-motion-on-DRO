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
load('FloquetEig12')
x0_DRO_0 = x0_DRO_M_3d;
x0_REL_period = -2/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';
x0_target = x0_DRO_0;

Flag_coordinate = 'LVLH';
% Flag_coordinate = 'VVLH';
% Flag_coordinate = 'VNC';

% isplot = 1;
isplot = 0;


t_total_all = [0.01:0.01:0.1,0.12:0.02:2]*para.T0;
% t_total_all = [0.01:0.01:0.1,0.12:0.02:1]*para.T0;
size_t = length(t_total_all);

% flag_str = 'Normal vector';
flag_str = '法向量';
r_normal_all = [0,0,1; 0,1,0; 1,0,0; 0,1,1];
size_r = size(r_normal_all,1);

% flag_str = 'Center';
% flag_str = '圆心';
% r_cent_all = [0,0,0; 0,0,1; 0,1,0; 1,0,0];
r_cent_all = [0,0,0; 1,0,0; 10,0,0; 100,0,0];
% r_cent_all = [0,0,0; 0,1,0; 0,10,0; 0,100,0];
% r_cent_all = [0,0,0; 0,0,1; 0,0,10; 0,0,100];
% r_cent_all = [0,0,0; 0.707,0.707,0; 7.07,7.07,0; 70.7,70.7,0];
% size_r = size(r_cent_all,1);

% flag_str = 'Radius';
% flag_str = '半径';
% r_radi_all = [1;10;100];
% size_r = size(r_radi_all,1);

patchflag = 1;

% flag_str = []; size_r = 1; patchflag = 2;
dv_total_m_all2 = zeros(size_t,size_r);
dv_total_xyz_all2 = zeros(size_t,3,size_r);
for ii_loop = 1:size_r
    %% 设定时间
    dv_total_m_all = zeros(size(t_total_all));
    dv_total_xyz_all = zeros(size_t,3);
    parfor kk_index = 1:size_t
        %% 选定受控方式
        
        if patchflag == 1 %% circumnavigation，绕飞
            n_samp = 10; % 每圈的点数
            dt = t_total_all(kk_index)/(n_samp-1); % 转移时间
    %         dt = 10/24/con.T_norma_day/(n_samp-1); % 转移时间
            % 空间圆方程
            if strcmp(flag_str,'法向量')
                r_normal = r_normal_all(ii_loop,:);
            else
                r_normal = [0,1,1]; % 法向量
            end
            r_normal = r_normal/norm(r_normal); % 单位法向量
            
            if strcmp(flag_str,'圆心')
                r_cent = r_cent_all(ii_loop,:);
            else
                r_cent = [0,0,0]; % 圆心，km
            end
            if strcmp(flag_str,'半径')
                radi = r_radi_all(ii_loop);
            else
                radi = 1; % 半径，km
            end
            
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
            
            n_samp = 6; % 变轨点数，若点数为2，则是两点间直接转移
    %         dt = dt_total/(n_samp-1);
            dt = t_total_all(kk_index)/(n_samp-1); % 转移时间

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
            plot3(r_chaser_km(:,1),r_chaser_km(:,2),r_chaser_km(:,3),'.'); axis equal
        else  %% 定点悬停
        %     dt = 24/24/con.T_norma_day; % 转移时间
        %     r_chaser_km = [3,3,3; ];
        end

        %% 变轨计算
        dt_all = dt*ones(size(r_chaser_km,1)-1,1); % 时间

        [x0_DRO_all,x0_REL_all_km,dv_all_km] = forcedRelMotion(x0_target,r_chaser_km,dt_all,Flag_coordinate,isplot,0);
        dvnorm_all_km = sqrt(sum(dv_all_km.^2,2));
        dv_total_m = sum(dvnorm_all_km)*1e3;

        dv_total_m_all(kk_index) = dv_total_m;
        dv_total_xyz_all(kk_index,:) = sum(abs(dv_all_km),1)*1e3;
        disp([kk_index, size_t])
    end
    dv_total_m_all2(:,ii_loop) = dv_total_m_all;
    dv_total_xyz_all2(:,:,ii_loop) = dv_total_xyz_all;
end
% semilogy(t_total_all*con.T_norma_day, dv_total_xyz_all2(:,:,1))

%% 画图
fig = figure(1);
Color_all = get(gca,'colororder');
p1 = semilogy(t_total_all*con.T_norma_day,dv_total_m_all2,'LineWidth',2);
for ii = 1:length(p1)
    p1(ii).Color = Color_all(ii,:);
end

% xlabel('Total time [day]'); ylabel('\Deltav [m/s]');
% title('Total \Deltav (CR3BP)');
xlabel('总时间 [day]'); ylabel('\Delta{\itv} [m/s]');
% title('总\Deltav');
grid on; grid minor
box on
view(0,90)
xlim([min([0,t_total_all*con.T_norma_day]),max(t_total_all*con.T_norma_day)])
% ylim([0.03,10^1])
% ylim([0.03,2*10^2])
fig.Renderer = 'painters';

str = strings(1,size_r);
for jj_loop = 1:size_r
    if strcmp(flag_str,'法向量')
        str(jj_loop) = ['',flag_str,'',' = [',num2str(r_normal_all(jj_loop,1)),',',...
            num2str(r_normal_all(jj_loop,2)),',',num2str(r_normal_all(jj_loop,3)),']^T'];
    elseif strcmp(flag_str,'圆心')
        str(jj_loop) = ['',flag_str,'',' = [',num2str(r_cent_all(jj_loop,1)),',',...
            num2str(r_cent_all(jj_loop,2)),',',num2str(r_cent_all(jj_loop,3)),']^T'];
    elseif strcmp(flag_str,'半径')
        str(jj_loop) = ['',flag_str,'',' = ',num2str(r_radi_all(jj_loop))];
    end
end
legend(str,'Location','northeast')

set(gcf,'Color',[255,255,255]/255);
% export_fig Impulse.png -r600