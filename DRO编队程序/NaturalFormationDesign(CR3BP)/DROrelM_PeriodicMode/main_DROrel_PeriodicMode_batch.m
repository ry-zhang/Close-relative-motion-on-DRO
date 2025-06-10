%% 三体轨道中的相对运动
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../../subF_eom(CR3BP)')

format longg
format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

%% 常数与变量
flag = '4-2';
load(['generalSolFFT_all',flag,'.mat'])
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
plot_flag = 1;
num = length(para);
xy_max_all = zeros(num,2);
time_interval = zeros(num,2);
angle_all = zeros(num,1);
% for i_index = 1:num
% for i_index = 11
for i_index = [1,100,200,300]
% for i_index = 300
    %% 周期相对运动
    dt = para(i_index).T0; % 积分时间
    length_t = 200;
    t_sample = linspace(0,dt,length_t);
    theta_sample = t_sample*2*pi/para(i_index).T0;
    % t_sample_day = t_sample*con.T_norma_day;
    relMot_L_linear0 = real(iDFTmatrix_theta(coe(i_index).N,0) * coe(i_index).c1_e3hat)';
    k0 = 1/(relMot_L_linear0(2)*con.r_norma); % 约束绕飞编队尺度
    c1_e3hat_k0 = k0*coe(i_index).c1_e3hat;
    relMot_L_linear = real(iDFTmatrix_theta(coe(i_index).N,theta_sample) * c1_e3hat_k0)';

    %% 画绝对运动轨道、相对运动在L及M坐标系中的轨道，及交叉点对应的位置
    %% 计算self-intersection points
    % 由于selfintersect函数只能计算小规模的向量，因此叠加做两次交叉点计算
    % 第一步，粗算，初步计算出交叉点的位置
    [~,~,intersectionPos1] = selfintersect(relMot_L_linear(1,:),relMot_L_linear(2,:));
    % 第二步，细算，在第一步得到的交叉点附近重复迭代计算
    t_temple1 = t_sample;
    t_temple2 = t_sample;
    relMot_L_temp1 = relMot_L_linear;
    intersectionPos2 = [intersectionPos1(1),intersectionPos1(2)+length_t];
    try
        while norm(relMot_L_temp1(1:3,intersectionPos2(1))-relMot_L_temp1(1:3,intersectionPos2(1)+1))>1e-10
            t_temple1 = linspace(t_temple1(intersectionPos2(1)),t_temple1(intersectionPos2(1)+1),length_t);
            t_temple2 = linspace(t_temple2(intersectionPos2(2)-length_t),t_temple2(intersectionPos2(2)-length_t+1),length_t);
            relMot_L_temp1 = real(iDFTmatrix_theta(coe(i_index).N,t_temple1*2*pi/para(i_index).T0) * c1_e3hat_k0)';
            relMot_L_temp2 = real(iDFTmatrix_theta(coe(i_index).N,t_temple2*2*pi/para(i_index).T0) * c1_e3hat_k0)';
            [x0,y0,intersectionPos2] = selfintersect([relMot_L_temp1(1,:),relMot_L_temp2(1,:)],[relMot_L_temp1(2,:),relMot_L_temp2(2,:)]);
        end
    catch
        intersectionPos2 = [intersectionPos1(1),intersectionPos1(2)+length_t];
    end

    % 得到两次交叉点对应的时间
    t_intersec(1) = t_temple1(intersectionPos2(1));
    t_intersec(2) = t_temple2(intersectionPos2(2)-length_t);
    rel_motion_L_linear_intersec = real(iDFTmatrix_theta(coe(i_index).N,t_intersec*2*pi/para(i_index).T0) * c1_e3hat_k0)';

    %% 将数据单位其转化为km
    rel_motion_L_linear_km = relMot_L_linear*con.r_norma;
    rel_motion_L_linear_intersec_km = rel_motion_L_linear_intersec*con.r_norma;

    %% 画图
    if plot_flag == 1
        figure(1)
        % x_ratio = 3.2; % x轴显示与图形中点的比例
        % y2x_ratio = 0.5; % 画图时y轴显示与x轴显示的比例
        x_ratio = 2; y2x_ratio = 1.2; hold off
        % plot(rel_motion_L_linear_km(1,1),rel_motion_L_linear_km(2,1),'g^'); hold on
%         plot(rel_motion_L_linear_intersec_km(1,:),rel_motion_L_linear_intersec_km(2,:),'k*'); hold on
        if max(rel_motion_L_linear_km(1,:)-1) <1e-4
            plot(rel_motion_L_linear_km(1,:),rel_motion_L_linear_km(2,:),'Color',[1, 114, 189]/255,'LineWidth',3); hold off
        else
            plot(rel_motion_L_linear_km(1,:),rel_motion_L_linear_km(2,:),'Color',[1, 114, 189]/255,'LineWidth',1.5); hold off
        end
        xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]');
        % legend('InitialPos','IntersecPos')
        axis equal; set(gca,'FontSize',17)
        grid on; grid minor
        x_max = max(rel_motion_L_linear_km(1,:));
        x_min = min(rel_motion_L_linear_km(1,:));
        x_middle = (x_max+x_min)/2;
        x_diff = x_max - x_min;
        y_middle = (max(rel_motion_L_linear_km(2,:))+min(rel_motion_L_linear_km(2,:)))/2;
        xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
        ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 
        xlim([-0.6,0.6]) 
%         xlim([-0.8,0.8]) 
        ylim([0,1.5]) 
        % title('周期相对运动 (坐标系L)')
        % legend('初始时刻','交叉点')
%         title('周期模态 (坐标系\itL\rm)')
%         title('平面周期模态')
        title(['\itT_{\rm0}\rm=',num2str(para(i_index).T0*con.T_norma_day,'%.2f'),'天'])
        % legend('Initial position','Intersection position')
%         legend('自相交点')
        set(gcf,'Color',[255,255,255]/255);
        exportgraphics(gcf,'PeMode.jpg','Resolution',600)
        pause(0.1)
    end
%     xy_ave = mean(abs(rel_motion_L_linear_km(1:2,:)),2);
%     xy_max_all(i_index,:) = max(abs(rel_motion_L_linear_km(1:2,:)-xy_ave),[],2)';
    xy_max_all(i_index,:) = max(rel_motion_L_linear_km(1:2,:),[],2)'-min(rel_motion_L_linear_km(1:2,:),[],2)';
    time_interval(i_index,:) = [diff(t_intersec),para(i_index).T0-diff(t_intersec)];
    angle_all(i_index) = 2*max(atan2d(rel_motion_L_linear_km(2,:),rel_motion_L_linear_km(1,:))-90);
%     计算视角，量化尺度（x方向、y方向宽度），内外圈时间比例？
end

%% 分析不同DRO的几何特性
figure(2)
subplot(1,2,1)
color_all = get(gca,'colororder');
plot([para.T0]*con.T_norma_day,xy_max_all,'LineWidth',1.5);
xlabel('\itT_{\rmDRO} \rm[day]'); ylabel('振幅 \rm[km]');
set(gca,'FontSize',15); xlim([0,27])
grid on; grid minor
legend('\itx_{\rmamp}','\ity_{\rmamp}','Location','north')
% xlim([-0.8,0.8]) 
% ylim([0,1.5]) 
subplot(1,2,2)
plot([para.T0]*con.T_norma_day,angle_all,'Color',color_all(3,:),'LineWidth',1.5);
xlabel('\itT_{\rmDRO} \rm[day]'); ylabel('\beta \rm[deg]');
set(gca,'FontSize',15); xlim([0,27])
grid on; grid minor


% set(gcf,'Color',[255,255,255]/255);
exportgraphics(gcf,'PeMode_all.jpg','Resolution',600)
% 存成png



