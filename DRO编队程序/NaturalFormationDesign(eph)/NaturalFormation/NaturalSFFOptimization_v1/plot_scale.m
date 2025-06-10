clear

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
%% 初始化
load('QDROforma1period.mat')
% load('DROforma1month.mat')
% load('DROforma2month.mat')

aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC  = t0UTC; % 初始历元
aux = initialize(aux); % 初始化

%% 计算优化结果（1倍放缩时的真值）
[~,num] = min(MaxDist_all);
x0_TC_optimize = x0_TC_optimize_all(num,:);
x0_TC_VVLH_chaser = [x0_TC_optimize(1),0,x0_TC_optimize(2),x0_TC_optimize(3)*1e-6,0,x0_TC_optimize(4)*1e-6];
% x0_MCR_target = xx_MCR_target(1,:);
[xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);
x0_MCR_chaser = T_TCO2TCR_eph(x0_TC_VVLH_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);
% 计算相对运动优化结果
rvTC_MCR = xx_MCR_chaser-xx_MCR_target;
rvTC_VVLH_ref = T_TCR2TCO_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'VVLH');


%% 计算不同放缩倍数时的末端误差
scale_all = [-2,-1,-0.5,0.5,1,2];
size_scale = length(scale_all);
error_all1 = zeros(size_scale,1);
error_all2 = zeros(size_scale,1);
figure(5)
for jj_index = 1:size_scale
% for jj_index = 1
    scale = scale_all(jj_index);
    x0_TC_optimize = scale*x0_TC_optimize_all(num,:);
    x0_TC_VVLH_chaser = [x0_TC_optimize(1),0,x0_TC_optimize(2),x0_TC_optimize(3)*1e-6,0,x0_TC_optimize(4)*1e-6];
    x0_MCR_chaser = T_TCO2TCR_eph(x0_TC_VVLH_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
    % 副星星历积分
    xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);
    % 计算相对运动
    rvTC_MCR = xx_MCR_chaser-xx_MCR_target;
    rv_TCR = T_TCR2TCO_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'LVLH');
    
    % 计算误差
    error_all1(jj_index) = norm(rv_TCR(end,1:3)-scale*rvTC_VVLH_ref(end,1:3));
    error_all2(jj_index) = sum(sqrt(sum((rv_TCR(:,1:3)-scale*rvTC_VVLH_ref(:,1:3)).^2,2)))/size(rv_TCR,1);
    
    % 画图
    plot(rv_TCR(:,1), rv_TCR(:,2),'LineWidth',2); hold on;
%     plot(rv_TCR(1,1), rv_TCR(1,2),'g^');
%     plot(rv_TCR(end,1), rv_TCR(end,2),'rv');
end

%% 图片设置
[~,num_min] = min(rv_TCR(:,1));
[~,num_max] = max(rv_TCR(:,1));
plot(1.5*[rv_TCR(num_min,1),-rv_TCR(num_min,1)],1.5*[rv_TCR(num_min,2),-rv_TCR(num_min,2)],'k--')
plot(1.5*[rv_TCR(num_max,1),-rv_TCR(num_max,1)],1.5*[rv_TCR(num_max,2),-rv_TCR(num_max,2)],'k--')
box on; grid on; grid minor; 
hold off; 
axis equal; xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]')
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
% title('Relative Orbit with different scales (LVLH)')
title('伴飞轨迹的线性放缩 (坐标系\itL)')
% title('Relative Orbit (LVLH, DRO1)')
set(gcf,'Renderer','painters')
str = strings(1,length(scale_all));
for ii = 1:length(scale_all)
    str(ii) = ['\itn \rm= ',num2str(scale_all(ii))];
end
% legend(str,'Location','north')
legend(str,'Location','east')