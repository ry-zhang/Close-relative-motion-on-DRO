%% 脉冲辅助下的DRO长期绕飞编队轨道设计与分析
% 2021-8-30
% by Yang Chihang
% email: ychhtl@foxmail.com

close all
clear
format longg
format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 13);%坐标轴
set(0,'defaultTextFontSize', 13);%文字
set(0,'defaultLineLineWidth',1.5)

addpath('../../subF_eom(CR3BP)')
load('generalSolFFT_12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

length_t = 201;

%% 在平面内拟周期轨道上寻找一段轨道，具有合适的距原点最小距离和首末端距离
num_ii = 16*10+1;
theta1_0_all = linspace(0,2*pi,num_ii); isplot = 0;
% theta1_0_all = linspace(0,0.8*pi,5); isplot = 1;
% theta1_0_all = [0.625,1.625]*pi; isplot = 1;
dt1 = para.T0;
t_sample2 = linspace(0,dt1,length_t);
e1_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e1hat)';
e2_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e2hat)';
dist_min = zeros(size(theta1_0_all));
dist_r0f = zeros(size(theta1_0_all));

for ii_loop = 1:length(theta1_0_all)
    
    
    theta1_0_temp = mod(theta1_0_all(ii_loop)+t_sample2*2*pi/para.T1,2*pi);
    rel_motion_temp = cos(theta1_0_temp).*e1_hat_refit + sin(theta1_0_temp).*e2_hat_refit;
    
    dist_min(ii_loop) = min(sqrt(sum(rel_motion_temp(1:2,:).^2)));
    dist_r0f(ii_loop) = min(sqrt(sum(diff(rel_motion_temp(1:2,[1,end])').^2)));
    if isplot
        figure(1)
        subplot(1,5,rem(ii_loop-1,5)+1)
        plot(rel_motion_temp(1,:),rel_motion_temp(2,:),'LineWidth',1.5); hold on
        plot(rel_motion_temp(1,1),rel_motion_temp(2,1),'g^')
        plot(rel_motion_temp(1,end),rel_motion_temp(2,end),'rv'); 
        plot(0,0,'ks');
        hold off
        axis equal; grid on
        xlim([-0.401,0.401]); ylim([-0.501,0.501])
        set(gca,'FontSize',13)
        xlabel('\itx_L/k_{\rm1}'); ylabel('\ity_L/k_{\rm1}');
        title(['\theta_{1,0}=',num2str(theta1_0_temp(1)/pi),'\pi'])
        pause(0.3)
    end
end
% legend('Trajectory','Initial Position','Final Position','Origin',...
%     'Location','southoutside', 'Orientation','horizontal')
legend('轨迹','初始位置','末端位置','主星',...
    'Location','southoutside', 'Orientation','horizontal')
% set(gcf,'Color',[255,255,255]/255);
% export_fig PBoundTrajS.png -r600
if isplot == 0
    figure(2)
    % subplot(1,2,1)
    plot(theta1_0_all,dist_min,'LineWidth',1.5);
    % set(gca,'FontSize',13); xlabel('\theta_{1,0}'); ylabel('minimum distance');
    % xlim([0,2*pi]); xticks([0,pi/2,pi,3*pi/2,2*pi]);
    % xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
    % subplot(1,2,2)
    hold on
    plot(theta1_0_all,dist_r0f,'LineWidth',1.5);
    set(gca,'FontSize',13); xlabel('\theta_{1,0} [rad]'); ylabel('distance');
    xlim([0,2*pi]); xticks([0,pi/2,pi,3*pi/2,2*pi]);
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
    ylim([0,1]); grid on; grid minor
    hold off
    le = legend('$\min \left\| {\rho }\left( t \right) \right\|$',...
        '$\left\| {\rho }\left( 0 \right) -{\rho }\left( T \right)\right\|$','Location','north');
    le.Interpreter = 'latex';
end
% set(gcf,'Color',[255,255,255]/255);
% export_fig PNaturalTraj.png -r600
%% 引入周期轨道(以及法向拟周期)，进一步增大最小距离

num_ii = 1000;
% theta1_0_all = [0.625,0.625,0.625,0.625,1.625,1.625,1.625,1.625]*pi;
% k0_all = [-0.1,-0.2,-0.3,-0.4,0.1,0.2,0.3,0.4];
k0_all = linspace(-1,0,num_ii);
theta1_0_all = 0.625*pi*ones(size(k0_all));
% k0_all = linspace(0,1,num_ii);
% theta1_0_all = 1.625*pi*ones(size(k0_all));
% theta2_0_all = [0.25,0.875,1.25,1.75]*pi;
dt1 = 1*para.T0;
t_sample2 = linspace(0,dt1,length_t);

e3_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e3hat)';
e1_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e1hat)';
e2_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e2hat)';
e5_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e5hat)';
e6_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e6hat)';

dist_min = zeros(size(theta1_0_all));
dist_r0f = zeros(size(theta1_0_all));


isplot = 0;
% k0_all = [0,-0.197,-0.743,-1];
% theta1_0_all = [0.625,0.625,0.625,0.625]*pi;
% isplot = 1;
for ii_loop = 1:length(theta1_0_all)

    theta1_0_temp = theta1_0_all(ii_loop);    
    theta1_temp = mod(theta1_0_temp+t_sample2*2*pi/para.T1,2*pi);
    rv_rel_planar = cos(theta1_temp).*e1_hat_refit + sin(theta1_temp).*e2_hat_refit;

    % 引入法向拟周期
%     theta2_0_temp = theta2_0_all(ii_loop);
%     theta2_temp = mod(theta2_0_temp+t_sample2*2*pi/para.T2,2*pi);
%     rv_rel_vertical = cos(theta2_temp).*e5_hat_refit + sin(theta2_temp).*e6_hat_refit;

    k0 = k0_all(ii_loop);
    rv_rel_periodic = k0*e3_hat_refit;
    
    rel_motion_t0f = rv_rel_planar + rv_rel_periodic;
    
    dist_min(ii_loop) = min(sqrt(sum(rel_motion_t0f(1:2,:).^2)));
    dist_r0f(ii_loop) = min(sqrt(sum(diff(rel_motion_t0f(1:2,[1,end])').^2)));
    
    if isplot
        % plot
        figure(5) 
        subplot(1,4,rem(ii_loop-1,4)+1)
%         plot3(rel_motion_t0f(1,:),rel_motion_t0f(2,:),rel_motion_t0f(3,:),'LineWidth',1.5); hold on
%         plot3(rel_motion_t0f(1,1),rel_motion_t0f(2,1),rel_motion_t0f(3,1),'g^')
%         plot3(rel_motion_t0f(1,end),rel_motion_t0f(2,end),rel_motion_t0f(3,end),'rv'); 
%         plot3(0,0,0,'ks');
        plot(rel_motion_t0f(1,:),rel_motion_t0f(2,:),'LineWidth',1.5); hold on
        plot(rel_motion_t0f(1,1),rel_motion_t0f(2,1),'g^')
        plot(rel_motion_t0f(1,end),rel_motion_t0f(2,end),'rv'); 
        plot(0,0,'ks'); hold off
%         view([0,90])
        
        grid on; axis equal; box on
        set(gcf,'Renderer','painters')
        xlim([-0.401,0.401]); ylim([-0.801,0.501]); zlim([-0.801,0.801]);
        set(gca,'FontSize',13)
        xlabel('\itx_L'); ylabel('\ity_L');
    %     title(['\theta_{1,0}=',num2str(theta1_0_temp/pi),'\pi, ','\theta_2(0)=',num2str(theta2_0_temp/pi),'\pi'])
%         title(['\theta_{1,0}=',num2str(theta1_0_temp/pi),'\pi,','k_0=',num2str(k0)])
        title(['\itk_{\rm0}\rm=',num2str(k0)])
        pause(0.3)
    end
end
legend('Trajectory','Initial Position','Final Position','Origin',...
    'Location','southoutside', 'Orientation','horizontal')

if isplot == 0
    figure(6)
    plot(k0_all,dist_min,'LineWidth',1.5); hold on
    plot(k0_all,dist_r0f,'LineWidth',1.5);
    set(gca,'FontSize',13); xlabel('\it{k}_{\rm0}'); ylabel('distance');
    ylim([0,0.3]); 
    grid on; grid minor
    hold off
    le = legend('$\min \left\| {\rho }\left( t \right) \right\|$',...
        '$\left\| {\rho }\left( 0 \right) -{\rho }\left( T \right)\right\|$','Location','north');
    le.Interpreter = 'latex';
%     legend('Minimum distance between spacecraft','Distance between initial and final positions','Location','north')
%     title(['\theta_{1,0}=',num2str(theta1_0_temp/pi),'\pi'])
%     set(gcf,'Color',[255,255,255]/255);
%     export_fig PPeDistance.png -r600
end

% export_fig PPeNaturalTraj.png -r600

%% 法向拟周期轨道一个周期间始末点距离分析
num_ii = 400;
theta2_0_all = linspace(0,2*pi,num_ii);
dt1 = 1*para.T0;
t_sample6 = linspace(0,dt1,length_t);
e5_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample6*2*pi/para.T0) * coe.c1_e5hat)';
e6_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample6*2*pi/para.T0) * coe.c1_e6hat)';
dist_r0f = zeros(size(theta2_0_all));


isplot = 1; 
theta2_0_all = [0.266,1.266]*pi;
% isplot = 0;
% theta2_0_all = linspace(0,pi/2,5);
for ii_loop = 1:length(theta2_0_all)
    
    theta2_temp = mod(theta2_0_all(ii_loop)+t_sample6*2*pi/para.T2,2*pi);
    rel_motion_temp = cos(theta2_temp).*e5_hat_refit + sin(theta2_temp).*e6_hat_refit;
    dist_r0f(ii_loop) = min(sqrt(sum(diff(rel_motion_temp(3,[1,end])').^2)));
    
    % plot
    if isplot
        figure(3);
        subplot(1,2,rem(ii_loop-1,2)+1)
        plot(t_sample6/para.T0,rel_motion_temp(3,:),'LineWidth',1.5)
        grid on
        xlim([0,dt1/para.T0]);
        ylim([-0.801,0.801]);
        set(gca,'FontSize',15)
        xlabel('{\itt} [{\itT}_0]'); ylabel('\itz_L');
        title(['\theta_{2,0}=',num2str(theta2_temp(1)/pi),'\pi'])
        pause(0.3)
    end
end
% set(gcf,'Color',[255,255,255]/255);
% export_fig NNaturalTrajectory.png -r600

if isplot == 0
    figure(4)
    plot(theta2_0_all,dist_r0f,'LineWidth',1.5);
    set(gca,'FontSize',15); xlabel('\theta_{2,0} [rad]'); ylabel('始末距离');
    xlim([0,2*pi]); xticks([0,pi/2,pi,3*pi/2,2*pi]);
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
    ylim([0,0.8]); grid on; grid minor
    % set(gcf,'Color',[255,255,255]/255);
    % export_fig NDistance.png -r600
end

%% 首末点变轨速度脉冲分析
sol = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 para.T0], para.x0_DRO, opts);

size_dt = 300;
dt_all = linspace(0.001,0.6,size_dt)*para.T0;

% size_dt = 4;
% dt_all = linspace(0.54,0.55,size_dt)*para.T0;
% dt_all = linspace(0.5461,0.5462,size_dt)*para.T0;

dv_total_m_all = zeros(1,size_dt);

theta1_0 = 0.625*pi; theta2_0 = 1.266*pi;
k0 = -0.197; k1 = 1; 
% k1 = ratio/con.r_norma; k0 = -0.197*k1;
k2 = 0.5; 
% k2 = 0; 

% isplot = 0;
% parfor jj_loop = 1:size_dt
% for jj_loop = 1
isplot = 1;
% dt_all = [0.001,0.01,0.1,0.2]*para.T0;
% dt_all = [0.3,0.4,0.5,0.6]*para.T0;
dt_all = [0.3,0.4,0.55,0.58]*para.T0;

for jj_loop = 1:4

    dt = dt_all(jj_loop);
    t0f = [para.T0-dt/2,dt/2];
    
    rel_motion_t0f = generalSol_relMotion(t0f,k0,k1,k2,theta1_0,theta2_0,para,coe);
    
    rel_motion_halfT = generalSol_relMotion(para.T0/2,k0,k1,k2,theta1_0,theta2_0,para,coe);
    ratio = 1/abs(rel_motion_halfT(2,1));
    rel_motion_t0f_km = rel_motion_t0f;
    rel_motion_t0f_km(1:3,:) = ratio*rel_motion_t0f(1:3,:);
    rel_motion_t0f_km(4:6,:) = ratio*rel_motion_t0f(4:6,:)*con.v_norma/con.r_norma;

    x0_target = deval(sol,t0f(1))'; 
    [x0_DRO_all,x0_REL_all_km,~,r_REL_traj_all_km] = forcedRelMotion(x0_target,rel_motion_t0f_km(1:3,:)',dt,'LVLH',0);
    dv_all_km = [x0_REL_all_km(1:2,4:6)-rel_motion_t0f_km(4:6,1:2)'];
    dvnorm_all_km = sqrt(sum(dv_all_km.^2,2));
    dv_total_m_all(jj_loop) = sum(dvnorm_all_km)*1e3;
    
    if isplot
        % 画图
        figure(7)
        subplot(1,4,rem(jj_loop-1,4)+1)
        t_sample3 = linspace(t0f(1),t0f(2),length_t);
        rel_motion_temp = generalSol_relMotion(t_sample3,k0,k1,k2,theta1_0,theta2_0,para,coe);
        rel_motion_temp_km = ratio*rel_motion_temp(1:3,:);
        p1 = plot3(rel_motion_temp_km(1,:),rel_motion_temp_km(2,:),rel_motion_temp_km(3,:),'LineWidth',1.5); hold on
        p2 = plot3(r_REL_traj_all_km(:,1),r_REL_traj_all_km(:,2),r_REL_traj_all_km(:,3),'LineWidth',1.5);
        p3 = plot3(rel_motion_t0f_km(1,1),rel_motion_t0f_km(2,1),rel_motion_t0f_km(3,1),'g^');
        p4 = plot3(rel_motion_t0f_km(1,2),rel_motion_t0f_km(2,2),rel_motion_t0f_km(3,2),'rv');
        p5 = plot3(0,0,0,'ks'); hold off
%         view([0,90])
%         view(0,0)
        view(-126,20)
%         view(-55,30)
        grid on; axis equal; box on
        set(gcf,'Renderer','painters')
        xlim([-0.801,0.801]); ylim([-1.201,1.201]); 
%         zlim([-0.5,0.4]);
        set(gca,'FontSize',16)
        label('\itx_L \rm[km]','\ity_L \rm[km]','\itz_L \rm[km]')
        title(['\Delta\itt\rm=',num2str(dt/para.T0),'{\itT}_0'])
        pause(0.3)
    end
end
% legend('Natural Trajectory','Transfer Trajectory','Initial Position','Final Position','Origin',...
%     'Location','southoutside', 'Orientation','vertical','NumColumns',3)
legend('自然轨迹','转移轨迹','自然轨迹起始位置','自然轨迹末端位置','主星',...
     'Location','southoutside', 'Orientation','vertical','NumColumns',3)
% set(gcf,'Color',[255,255,255]/255);
% export_fig 3DFlyAroundTraj2.png -r600

if isplot == 0
    figure(8)
    semilogy(dt_all/para.T0,dv_total_m_all,'LineWidth',1.5)
    ylim([1e-3,1e0])
    set(gca,'FontSize',15); xlabel('\Delta{\itt} [{\itT}]'); ylabel('\Delta\itv \rm[m/s]');
    grid on; grid minor
    set(gcf,'Color',[255,255,255]/255);
    export_fig TransferImp3D.png -r600
end
