close all
clear
addpath('../../subF_eom(CR3BP)')
load('generalSolFFT_12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 15);%坐标轴
set(0,'defaultTextFontSize', 15);%文字
set(0,'defaultLineLineWidth',1.5)

%% 计算ECI下的相对运动
dt = 2*para.T0; % 积分时间                  
length_t = 2000;
delta_t = 0.2*para.T0;
t_sample = linspace(delta_t/2,para.T0-delta_t/2,length_t);
t_sample_day = t_sample*con.T_norma_day;
ratio = 1;

% 自然轨迹 绝对运动
sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 dt], para.x0_DRO, opts);
x_chief_MCR = deval(sol_DRO,t_sample);
% 将月球为中心天体的旋转坐标系转化到地球为中心天体的旋转坐标系，绝对运动转换
x_chief_ECR = x_chief_MCR; x_chief_ECR(1:3,:) = x_chief_MCR(1:3,:)-[-1,0,0]';
x_chief_ECI = synodic2inertial(x_chief_ECR,t_sample);

k1_all = [1,0.5,0.2,-0.2,-0.5,-1]*1e-3;
k2_all = [0.2,0.5,0.2,-0.2,-0.3,0.3]*1e-3;
for ii_loop = 1:length(k1_all)
    % 自然轨迹 相对运动
    theta1_0 = 0.625*pi; theta2_0 = 1.266*pi;
    k1 = k1_all(ii_loop); 
    k2 = k2_all(ii_loop);
%     k2 = 0.1; 
    k0 = -0.197*k1;
    x_rel_L_natural = generalSol_relMotion(t_sample,k0,k1,k2,theta1_0,theta2_0,para,coe);
    x_rel_MCR_natural = T_TCO2TCR_CR3BP(x_rel_L_natural',x_chief_MCR','LVLH',con.mu)';

    % 受控轨迹 绝对运动/相对运动
    x0_chief_MCR_transfer = deval(sol_DRO,t_sample(end))'; 
    [x0f_chief_MCR_transfer,x0f_rel_L_transfer_km] = forcedRelMotion(x0_chief_MCR_transfer,x_rel_L_natural(1:3,[end,1])',delta_t,'LVLH',0,0);
    x0_rel_L_transfer = [x0f_rel_L_transfer_km(1,1:3),x0f_rel_L_transfer_km(1,4:6)*con.r_norma/con.v_norma];
    % t_sample2 = linspace(para.T0-delta_t/2,para.T0+delta_t/2,length_t/2);
    t_sample2 = linspace(0,delta_t,length_t/2);
    sol_transfer = ode113(@(t,x)eom_rel3b(t,x,con.mu),t_sample2([1,end]), [x0f_chief_MCR_transfer(1,:), x0_rel_L_transfer], opts);
    sol_sample2 = deval(sol_transfer,t_sample2);
    x_chief_MCR_transfer = sol_sample2(1:6,:);
    x_rel_L_transfer = sol_sample2(7:12,:);
    x_rel_MCR_transfer = T_TCO2TCR_CR3BP(x_rel_L_transfer',x_chief_MCR_transfer','LVLH',con.mu)';

    % 将月球为中心天体的旋转坐标系转化到地球为中心天体的旋转坐标系，绝对运动转换
    x_chief_ECR_transfer= x_chief_MCR_transfer; x_chief_ECR_transfer(1:3,:) = x_chief_MCR_transfer(1:3,:)-[-1,0,0]';
    x_chief_ECI_transfer = synodic2inertial(x_chief_ECR_transfer,t_sample2);

    % 放缩
    r_rel_L_natural = ratio*con.r_norma*x_rel_L_natural(1:3,:);
    r_rel_L_transfer = ratio*con.r_norma*x_rel_L_transfer(1:3,:);
    
    % 地球为中心天体的惯性坐标系 下的轨道坐标系
    % natural 轨道段
    x_rel_ECR_natural = x_rel_MCR_natural; % 旋转坐标系下的相对运动是相同的，与旋转坐标系质心无关
    x_rel_ECI_natural = synodic2inertial(x_rel_ECR_natural,t_sample);
    x_rel_ECIL_natural = T_TCI2TCO_E_CR3BP(x_rel_ECI_natural',x_chief_ECI',t_sample,'LVLH',con.mu)';
    % transfer 轨道段
    x_rel_ECR_transfer = x_rel_MCR_transfer; % 旋转坐标系下的相对运动是相同的，与旋转坐标系质心无关
    x_rel_ECI_transfer = synodic2inertial(x_rel_ECR_transfer,t_sample2);
    x_rel_ECIL_transfer = T_TCI2TCO_E_CR3BP(x_rel_ECI_transfer',x_chief_ECI_transfer',t_sample2,'LVLH',con.mu)';
    % 由于DRO绕月逆行，因此LVLH的z轴方向与MCR的z轴相反。
    
    r_rel_ECIL_natural_km = ratio*con.r_norma*x_rel_ECIL_natural(1:3,:);
    r_rel_ECIL_transfer_km = ratio*con.r_norma*x_rel_ECIL_transfer(1:3,:);
    
    % 画MCR LVLH中的编队
    f1 = figure(2);
    f1.Position = [1000,340,455,420];
    plot3(r_rel_L_natural(1,:),r_rel_L_natural(2,:),r_rel_L_natural(3,:),...
        'LineWidth',1.5,'Color',[0, 114, 189]/255);
    hold on
    plot3(r_rel_L_transfer(1,:),r_rel_L_transfer(2,:),r_rel_L_transfer(3,:),...
        'LineWidth',1.5,'Color',[217, 81, 22]/255);
    plot3(r_rel_L_natural(1,1),r_rel_L_natural(2,1),r_rel_L_natural(3,1),'g^','MarkerSize',3);
    plot3(r_rel_L_natural(1,end),r_rel_L_natural(2,end),r_rel_L_natural(3,end),'rv','MarkerSize',3);
    plot(0,0,'ks','MarkerSize',5)
    xlabel('\itx_L {\rm[km]}'); ylabel('\ity_L {\rm[km]}'); 
    hold on; axis equal; grid on; box on
    xlim([min(k1_all),max(k1_all)]*con.r_norma*0.45);
    ylim([min(k1_all),max(k1_all)]*con.r_norma*0.35);
    view([0,90])
    set(gcf,'Color',[255,255,255]/255);
%     export_fig ImpulseFormMCR.png -r600
    
    % 画ECI LVLH中的编队
    f2 = figure(3);
    f2.Position = [1500,340,455,420];
    plot3(r_rel_ECIL_natural_km(1,:),r_rel_ECIL_natural_km(2,:),r_rel_ECIL_natural_km(3,:),...
        'LineWidth',1.5,'Color',[0, 114, 189]/255);
    hold on
    plot3(r_rel_ECIL_transfer_km(1,:),r_rel_ECIL_transfer_km(2,:),r_rel_ECIL_transfer_km(3,:),...
        'LineWidth',1.5,'Color',[217, 81, 22]/255);
    plot3(r_rel_ECIL_natural_km(1,1),r_rel_ECIL_natural_km(2,1),r_rel_ECIL_natural_km(3,1),'g^','MarkerSize',3);
    plot3(r_rel_ECIL_natural_km(1,end),r_rel_ECIL_natural_km(2,end),r_rel_ECIL_natural_km(3,end),'rv','MarkerSize',3);
    plot(0,0,'ks','MarkerSize',5)
    xlabel('\itx_{LE} {\rm[km]}'); ylabel('\ity_{LE} {\rm[km]}'); zlabel('\itz_{LE} {\rm[km]}'); 
    hold on; axis equal; grid on; box on
    xlim([min(k1_all),max(k1_all)]*con.r_norma*0.35);
    ylim([min(k1_all),max(k1_all)]*con.r_norma*0.35);
%     view([0,90])
    view([61,27])
    set(gcf,'Color',[255,255,255]/255);
%     export_fig ImpulseFormECI.png -r600
    
end

figure(2)
% legend('自然轨道段','转移轨道段','自然轨道段起始点','自然轨道段末端点','主航天器',...
%      'Location','southoutside', 'Orientation','vertical','NumColumns',3)
% set(gcf,'Color',[255,255,255]/255);
export_fig ImpulseFormMCR.png -r600

figure(3)
% % legend('自然轨道段','转移轨道段','自然轨道段起始点','自然轨道段末端点','主航天器',...
% %      'Location','northeastoutside')
% set(gcf,'Color',[255,255,255]/255);
export_fig ImpulseFormECI.png -r600
