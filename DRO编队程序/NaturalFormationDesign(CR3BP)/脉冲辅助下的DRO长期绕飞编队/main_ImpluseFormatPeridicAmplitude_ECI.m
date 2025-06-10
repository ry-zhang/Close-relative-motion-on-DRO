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
ratio = 5000/con.r_norma;

% 自然轨迹 绝对运动
sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 dt], para.x0_DRO, opts);
x_chief_MCR = deval(sol_DRO,t_sample);
% 将月球为中心天体的旋转坐标系转化到地球为中心天体的旋转坐标系，绝对运动转换
x_chief_ECR = x_chief_MCR; x_chief_ECR(1:3,:) = x_chief_MCR(1:3,:)-[-1,0,0]';
x_chief_ECI = synodic2inertial(x_chief_ECR,t_sample);

% 预分配存储轨迹的矩阵zry
x_rel_L_force = [];%zry

% k0_all = [0,-0.197,-0.743,-1];
k0_all = [-0.197];
for ii_loop = 1:length(k0_all)
    % 自然轨迹 相对运动
    theta1_0 = 0.625*pi; theta2_0 = 1.266*pi;
    k1 = 1; 
    k2 = 0;
%     k2 = 0.1; 
    k0 = k0_all(ii_loop);
    x_rel_L_natural = generalSol_relMotion(t_sample,k0,k1,k2,theta1_0,theta2_0,para,coe);
    x_rel_MCR_natural = T_TCO2TCR_CR3BP(x_rel_L_natural',x_chief_MCR','LVLH',con.mu)';

    % 受控轨迹 绝对运动/相对运动
    x0_chief_MCR_transfer = deval(sol_DRO,t_sample(end))'; 
    [x0f_chief_MCR_transfer,x0f_rel_L_transfer_km] = forcedRelMotion(x0_chief_MCR_transfer,x_rel_L_natural(1:3,[end,1])',delta_t,'LVLH',1,0);
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

%% x_rel_MCR_force 将存储相位角（单位为度）以及相对运动轨迹的 x, y, z, vx, vy, vz 信息，并被保存为 x_rel_MCR_force.mat。
% 计算相位 zry
    phase_natural = linspace(0, 180, size(x_rel_L_natural, 2)); % 自然轨道段
    phase_transfer = linspace(180, 360, size(x_rel_L_transfer, 2)); % 转移轨道段

    % 组合轨迹数据
    x_rel_L_force = [x_rel_L_force; ...
        [phase_natural',  x_rel_L_natural']; ...
        [phase_transfer', x_rel_L_transfer']];

    %% 画MCR LVLH中的编队
    f1 = figure(1);
    f1.Position = [1000,340,400,300];
    %mcr  %zry
%     plot3(x_rel_MCR_natural(1,:),x_rel_MCR_natural(2,:),x_rel_MCR_natural(3,:),'LineWidth',1.5);
%     hold on
%     plot3(x_rel_MCR_transfer(1,:),x_rel_MCR_transfer(2,:),x_rel_MCR_transfer(3,:),'LineWidth',1.5);
%     plot3(x_rel_MCR_natural(1,1),x_rel_MCR_natural(2,1),x_rel_MCR_natural(3,1),'g^');
%     plot3(x_rel_MCR_natural(1,end),x_rel_MCR_natural(2,end),x_rel_MCR_natural(3,end),'rv');
    %lvlh %zry
    plot3(r_rel_L_natural(1,:),r_rel_L_natural(2,:),r_rel_L_natural(3,:),'LineWidth',1.5);
    hold on
    plot3(r_rel_L_transfer(1,:),r_rel_L_transfer(2,:),r_rel_L_transfer(3,:),'LineWidth',1.5);
    plot3(r_rel_L_natural(1,1),r_rel_L_natural(2,1),r_rel_L_natural(3,1),'g^');
    plot3(r_rel_L_natural(1,end),r_rel_L_natural(2,end),r_rel_L_natural(3,end),'rv');
    plot(0,0,'ks','MarkerSize',5)
    xlabel('\itx_L'); ylabel('\ity_L'); 
    hold off; axis equal; grid on; box on
%     title({'脉冲辅助绕飞编队','(MCR LVLH)'}); 
    xlim([-0.401,0.401]); ylim([-0.801,0.501]); zlim([-0.801,0.801]);
%     title({['\itk_{\rm0}\rm=',num2str(k0)],'(MCR LVLH)'})
    view([0,90])
    set(gcf,'Color',[255,255,255]/255);
%     export_fig ImpulseFormMCR.png -r600
    
    % 画ECI LVLH中的编队
    f2 = figure(2);
    f2.Position = [1500,340,400,300];
    plot3(r_rel_ECIL_natural_km(1,:),r_rel_ECIL_natural_km(2,:),r_rel_ECIL_natural_km(3,:),'LineWidth',1.5);
    hold on
    plot3(r_rel_ECIL_transfer_km(1,:),r_rel_ECIL_transfer_km(2,:),r_rel_ECIL_transfer_km(3,:),'LineWidth',1.5);
    plot3(r_rel_ECIL_natural_km(1,1),r_rel_ECIL_natural_km(2,1),r_rel_ECIL_natural_km(3,1),'g^');
    plot3(r_rel_ECIL_natural_km(1,end),r_rel_ECIL_natural_km(2,end),r_rel_ECIL_natural_km(3,end),'rv');
    plot(0,0,'ks','MarkerSize',5)
    xlabel('\itx_L_E'); ylabel('\ity_L_E'); 
    hold off; axis equal; grid on;box on
%     title({'脉冲辅助绕飞编队','(MCR LVLH)'}); 
    xlim([-0.401,0.401]); ylim([-0.801,0.501]); zlim([-0.801,0.801]);
%     title({['\itk_{\rm0}\rm=',num2str(k0)],'(ECI LVLH)'})
    view([0,90])
    set(gcf,'Color',[255,255,255]/255);
% %     export_fig ImpulseFormECI.png -r600
    
end
legend('自然轨道段','转移轨道段','自然轨道段起始点','自然轨道段末端点','主航天器',...
     'Location','southoutside', 'Orientation','vertical','NumColumns',3)

% 保存数据到MAT文件zry
save('x_rel_L_force.mat', 'x_rel_L_force');
% figure(1)
% legend('自然轨迹','转移轨迹','自然轨迹起始位置','自然轨迹末端位置','主星',...
%      'Location','northeastoutside')
% set(gcf,'Color',[255,255,255]/255);
% export_fig ImpulseFormMCR.png -r600

% figure(2)
% legend('自然轨迹','转移轨迹','自然轨迹起始位置','自然轨迹末端位置','主星',...
%      'Location','northeastoutside')
% set(gcf,'Color',[255,255,255]/255);
% export_fig ImpulseFormECI.png -r600
