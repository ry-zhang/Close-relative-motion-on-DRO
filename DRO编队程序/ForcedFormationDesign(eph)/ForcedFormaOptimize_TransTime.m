%% 星历下的强迫相对运动
% 2021-7-11
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../subF_eom(CR3BP)')
addpath('../subF_eom(eph)')

format longg
format compact
warning off

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'DefaultAxesFontsize',15);
set(0,'DefaultTextFontsize',15);
set(0,'defaultLineLineWidth',1.5)

load('FloquetEig12.mat')

%% 
load('DROephMultiRev2024.mat')

% 加载星历、设置初始历元
aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;

% 星历积分初值（旋转系）
aux.jd0 = jd0;
aux = initialize(aux); % 初始化
x0_MCR_target = T_ECJ2k2Rot(aux.jd0, x0j2k', [0,0,0],aux.C_Mat, 'MCEMR');
% 星历积分初值（地月旋转系）
% x0_MCR_target = [72687.2175909459 0 0 0 -0.540957166578942 0]; % 星历DRO初值

% load('data_forced.mat')
% aux.jd0  = data.jd; % 初始历元
% aux.t0UTC = [];
% x0_MCR_target = data.xMT_MCR; % 星历DRO初值


%% 

Flag_coordinate = 'LVLH';
% Flag_coordinate = 'VVLH';
% Flag_coordinate = 'VNC';

% isplot = 1;
isplot = 0;
isdisplay = 0;

% t_total_all = [0.01:0.01:0.1,0.2:0.05:2]*para.T0*con.T_norma;
t_total_all = [0.01:0.01:0.1,0.12:0.02:2]*para.T0*con.T_norma;
size_t = length(t_total_all);

% r_normal_all = [0,0,1];
flag_str = '法向量';
r_normal_all = [0,0,1; 0,1,0; 1,0,0; 0,1,1];
size_r = size(r_normal_all,1);

% flag_str = '圆心';
% r_cent_all = [0,0,0; 0,0,1; 0,1,0; 1,0,0];
r_cent_all = [0,0,0; 1,0,0; 10,0,0; 100,0,0];
% r_cent_all = [0,0,0; 0,1,0; 0,10,0; 0,100,0];
% r_cent_all = [0,0,0; 0,0,1; 0,0,10; 0,0,100];
% r_cent_all = [0,0,0; 0.707,0.707,0; 7.07,7.07,0; 70.7,70.7,0];
% size_r = size(r_cent_all,1);

% flag_str = '半径';
% r_radi_all = [1;10;100];
% size_r = size(r_radi_all,1);

dv_total_m_all2 = zeros(size_t,size_r);
dv_total_xyz_all2 = zeros(size_t,3,size_r);
for ii_loop = 1:size_r
    %% 设定时间
    dv_total_m_all = zeros(size(t_total_all));
    dv_total_xyz_all = zeros(size_t,3);
    parfor kk_index = 1:size_t
        %% 选定受控方式
        patchflag = 1;
        if patchflag == 1 %% circumnavigation，绕飞
            n_samp = 10; % 每圈的点数
            dt = t_total_all(kk_index)/(n_samp-1); % 转移时间
    %         dt = 10/24/para.T_norma_day/(n_samp-1); % 转移时间
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
            r0 = [1,0,0];
            v0 = [0,0,0];
            rf = [2,0,0];
            
            n_samp = 2; % 变轨点数，若点数为2，则是两点间直接转移
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
        %     dt = 24/24/para.T_norma_day; % 转移时间
        %     r_chaser_km = [3,3,3; ];
        end

    %% 星历积分
    dt_all = dt*ones(size(r_chaser_km,1)-1,1); % 时间

    tic;
    [x0_DRO_all,x0_REL_all,dv_all] = forcedRelMotion_eph(x0_MCR_target,r_chaser_km,dt_all,aux,Flag_coordinate,isplot,isdisplay);
    time = toc;
    dvnorm_all = sqrt(sum(dv_all.^2,2));
    dv_total_m = sum(dvnorm_all)*1e3;

    % 如果初始状态的速度不是0，则在初始时刻的变轨脉冲上减去v0
    if patchflag~=1
        dv_all(1,:) = dv_all(1,:)-v0;
    end

    % 计算变轨点的jd
    jd_all = zeros(n_samp,1);
    jd_all(1) = aux.jd0;
    for jj_loop = 2:n_samp
        jd_all(jj_loop) = aux.jd0+sum(dt_all(1:jj_loop-1))/86400;
    end

    % save ForcedTransfer jd_all x0_DRO_all x0_REL_all dv_all

        dv_total_m_all(kk_index) = dv_total_m;
        dv_total_xyz_all(kk_index,:) = sum(abs(dv_all),1)*1e3;
        disp([kk_index, size_t])
    end
    dv_total_m_all2(:,ii_loop) = dv_total_m_all;
    dv_total_xyz_all2(:,:,ii_loop) = dv_total_xyz_all;
end

warning on
%% 画图

fig = figure(1); hold on
Color_all = get(gca,'colororder');
p2 = semilogy(t_total_all/86400,dv_total_m_all2,'LineWidth',2,'LineStyle',':');
for ii = 1:length(p2)
    p2(ii).Color = Color_all(ii,:);
end
xlabel('总时间 [day]'); ylabel('\Delta{\itv} [m/s]');
% title('Total \Deltav (eph)');
% grid on; grid minor
box on
view(0,90)
% xlim([min(t_total_all/86400),max(t_total_all/86400)])
xlim([min([0,t_total_all/86400]),max(t_total_all/86400)])
% ylim([0.01,10^1])
fig.Renderer = 'painters';

str = strings(1,size_r);
for jj_loop = 1:size_r
    if strcmp(flag_str,'法向量')
        str(jj_loop) = ['[',num2str(r_normal_all(jj_loop,1)),',',...
            num2str(r_normal_all(jj_loop,2)),',',num2str(r_normal_all(jj_loop,3)),']^T','(CRTBP)'];
    elseif strcmp(flag_str,'圆心')
        str(jj_loop) = ['[',num2str(r_cent_all(jj_loop,1)),',',...
            num2str(r_cent_all(jj_loop,2)),',',num2str(r_cent_all(jj_loop,3)),']^T','(CRTBP)'];
    elseif strcmp(flag_str,'半径')
        str(jj_loop) = [num2str(r_radi_all(jj_loop)),'(CRTBP)'];
    end
end

for jj_loop = 1:size_r
    if strcmp(flag_str,'法向量')
        str(jj_loop+size_r) = ['[',num2str(r_normal_all(jj_loop,1)),',',...
            num2str(r_normal_all(jj_loop,2)),',',num2str(r_normal_all(jj_loop,3)),']^T','(星历)'];
    elseif strcmp(flag_str,'圆心')
        str(jj_loop+size_r) = ['[',num2str(r_cent_all(jj_loop,1)),',',...
            num2str(r_cent_all(jj_loop,2)),',',num2str(r_cent_all(jj_loop,3)),']^T','(星历)'];
    elseif strcmp(flag_str,'半径')
        str(jj_loop+size_r) = [num2str(r_radi_all(jj_loop)),'(星历)'];
    end
end
legend(str,'Location','northeast','NumColumns',2)

set(gcf,'Color',[255,255,255]/255);
export_fig Impulse.png -r600