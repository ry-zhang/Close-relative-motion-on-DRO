clear
error_all1 = [];
for ii_loop = 1:3
    %% 初始化
    if ii_loop == 1
        load('QDROforma1period.mat')
    elseif ii_loop == 2
        load('QDROforma2period.mat')
    else
        load('QDROforma4period.mat')
    end

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
    rv_TCR_ref = T_TCR2TCO_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'VVLH');


    %% 计算不同放缩倍数时的末端误差
    % scale_all = [-10:1:-2,-1:0.1:1,2:10];
    scale_all = [-400:3:-2,-1:0.1:1,3:400];
    size_scale = length(scale_all);
    error_temp1 = zeros(size_scale,1);
%     error_temp2 = zeros(size_scale,1);

    parfor jj_index = 1:size_scale
    % for jj_index = 1
        scale = scale_all(jj_index);
        x0_TC_optimize = scale*x0_TC_optimize_all(num,:);
        x0_TC_VVLH_chaser = [x0_TC_optimize(1),0,x0_TC_optimize(2),x0_TC_optimize(3)*1e-6,0,x0_TC_optimize(4)*1e-6];
        x0_MCR_chaser = T_TCO2TCR_eph(x0_TC_VVLH_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
        % 副星星历积分
        xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);
        % 计算相对运动
        rvTC_MCR = xx_MCR_chaser-xx_MCR_target;
        rv_TCR = T_TCR2TCO_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'VVLH');

        % 计算误差
        error_temp1(jj_index) = norm(rv_TCR(end,1:3)-scale*rv_TCR_ref(end,1:3));
%         error_temp2(jj_index) = sum(sqrt(sum((rv_TCR(:,1:3)-scale*rv_TCR_ref(:,1:3)).^2,2)))/size(rv_TCR,1);
    end
    error_all1 = [error_all1, error_temp1];
end

%% 画绝对误差
figure(15)
plot(scale_all,error_all1(:,1),'LineWidth',1.5); hold on
plot(scale_all,error_all1(:,2),'LineWidth',1.5); 
plot(scale_all,error_all1(:,3),'LineWidth',1.5); 
ylabel('\ite\rm(\itt_f\rm) [km]')

box on; grid on; grid minor; hold off
xlabel('\itnS \rm[km]'); 

set(gca,'FontSize',15,'fontname','times new roman'); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
legend('0.5 Month','1 Month','2 Month','Location','north')
legend('13.64 days','27.28 days','54.57 days','Location','north')

%% 画相对误差
figure(16)
plot(scale_all,error_all1(:,1)./abs(scale_all'),'LineWidth',1.5); hold on
plot(scale_all,error_all1(:,2)./abs(scale_all'),'LineWidth',1.5); 
plot(scale_all,error_all1(:,3)./abs(scale_all'),'LineWidth',1.5); 
ylabel('\ite\rm(\itt_f\rm)/\it/S\rm/\itn')

box on; grid on; grid minor; hold off
xlabel('\itnS \rm[km]'); 

set(gca,'FontSize',15,'fontname','times new roman'); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
legend('0.5 Month','1 Month','2 Month','Location','north')
legend('13.64 days','27.28 days','54.57 days','Location','north')