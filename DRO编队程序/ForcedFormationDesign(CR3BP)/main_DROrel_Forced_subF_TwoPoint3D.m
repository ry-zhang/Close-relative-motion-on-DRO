%% CR3BP下的强迫相对运动
% 2021-7-11
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../subF_eom(CR3BP)')

format longg
format compact
close all
clear

addpath('../subF_eom(CR3BP)')
load('generalSolFFT_12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 15);%坐标轴
set(0,'defaultTextFontSize', 15);%文字
set(0,'defaultLineLineWidth',1.5)

%% 常数与变量

Flag_coordinate = 'LVLH';
% Flag_coordinate = 'VVLH';
% Flag_coordinate = 'VNC';

isplot = 1;
% isplot = 0;


% 星历下呢？

sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 para.T0], para.x0_DRO, opts);

%% 选择patch pointd 的方式
% k00 = 8.624e-6; k20 = 0; theta2_00 = 0;
% k0f = -8.624e-6; k2f = 0; theta2_0f = 0;
k00 = 5e-5; k20 = 0; theta200 = 0.5;
k0f = -5e-5; k2f = 0; theta20f = 0.5;
k1 = 0; theta1_0 = 0.625*pi; 

t_sample = linspace(0,para.T0,200);
x0_parking = generalSol_relMotion(t_sample,k00,k1,k20,theta1_0,theta200,para,coe);
xf_parking = generalSol_relMotion(t_sample,k0f,k1,k2f,theta1_0,theta20f,para,coe);
r0_parking_km = x0_parking(1:3,:)'*con.r_norma;
rf_parking_km = xf_parking(1:3,:)'*con.r_norma;

t0 = 5/con.T_norma_day;
dt = 3/con.T_norma_day;
tf = t0+dt;

x0 = generalSol_relMotion(t0,k00,k1,k20,theta1_0,theta200,para,coe);
xf = generalSol_relMotion(tf,k0f,k1,k2f,theta1_0,theta20f,para,coe);
x0_DRO = deval(sol_DRO,t0);
rv_chaser_km = [x0(1:3)*con.r_norma,xf(1:3)*con.r_norma; 
    x0(4:6)*con.v_norma,xf(4:6)*con.v_norma]';
[x0_DRO_all,r0_REL_all_km,~,r_REL_traj_all_km] = forcedRelMotion(x0_DRO',rv_chaser_km,dt,Flag_coordinate,0,0);
size_r_traj = size(r_REL_traj_all_km,1);
r_norm_REL_traj_all_km = sqrt(sum(r_REL_traj_all_km(:,1:2).^2,2));
[r_min_d,r_min_index] = min(r_norm_REL_traj_all_km);
t_min = t0+dt/size_r_traj*r_min_index;


k2 = 10e-6;

theta200_all = linspace(0,2*pi,101);
theta20f_all = linspace(0,2*pi,99);
dv_total_m_all = zeros(101,99);
r_min_normal_all = dv_total_m_all;


for theta0_loop = 1:length(theta200_all)
    parfor thetaf_loop = 1:length(theta20f_all)
        isplot = 0;
        theta200 = theta200_all(theta0_loop);
        theta20f = theta20f_all(thetaf_loop);
        
% for theta0_loop = 1:length(theta200_all)
%     for thetaf_loop = 1:length(theta20f_all)
%         isplot = 1;
%         theta200 = pi;
%         theta20f = pi;
        
        
        x0 = generalSol_relMotion(t0,k00,0,k2,0,theta200,para,coe);
        xf = generalSol_relMotion(tf,k0f,0,k2,0,theta20f,para,coe);
        
        x0_DRO = deval(sol_DRO,t0);
        rv_chaser_km = [x0(1:3)*con.r_norma,xf(1:3)*con.r_norma; 
        x0(4:6)*con.v_norma,xf(4:6)*con.v_norma]';
        [x0_DRO_all,r0_REL_all_km,~,r_REL_traj_all_km] = forcedRelMotion(x0_DRO',rv_chaser_km,dt,Flag_coordinate,0,0);

        r_min_p = r_REL_traj_all_km(r_min_index,:);
        % 记录最短距离
        r_min_normal_all(theta0_loop,thetaf_loop) = r_min_p(3);

        dv_all = r0_REL_all_km(:,4:6) - rv_chaser_km(:,4:6);
        dv_all(2,:) = -dv_all(2,:);
        dvnorm_all_km = sqrt(sum(dv_all.^2,2));
        dv_total_m = sum(dvnorm_all_km)*1e3;
        dv_total_m_all(theta0_loop,thetaf_loop) = dv_total_m;
        
        if isplot == 1
            x0_parking = generalSol_relMotion(t_sample,k00,0,k2,0,theta200,para,coe);
            xf_parking = generalSol_relMotion(t_sample,k0f,0,k2,0,theta20f,para,coe);
            r0_parking_km = x0_parking(1:3,:)'*con.r_norma;
            rf_parking_km = xf_parking(1:3,:)'*con.r_norma;
            figure(1)
            colorTable = get(gca,'colororder');
            % hold on
            p1 = plot3(0,0,0,'ks');  hold on
            p2 = plot3(r0_REL_all_km(1,1),r0_REL_all_km(1,2),r0_REL_all_km(1,3),'g^');
            p3 = plot3(r0_REL_all_km(end,1),r0_REL_all_km(end,2),r0_REL_all_km(end,3),'rv'); 

            p4 = plot3(r_REL_traj_all_km(:,1),r_REL_traj_all_km(:,2),r_REL_traj_all_km(:,3),'r'); 
            % n_samp_ref = 180;
            % theta_all = linspace(0,2*pi*n_samp_ref/(n_samp_ref+1),n_samp_ref);
            % r_chaser_km_ref = r_cent + radi*cos(theta_all)'*r_a + radi*sin(theta_all)'*r_b;
            p5 = plot3(r0_parking_km(:,1),r0_parking_km(:,2),r0_parking_km(:,3),'Color',colorTable(1,:));
            p6 = plot3(rf_parking_km(:,1),rf_parking_km(:,2),rf_parking_km(:,3),'Color',colorTable(1,:)); axis equal
            p71 = plot3(r_min_p(1),r_min_p(2),r_min_p(3),'ok'); axis equal
            p72 = plot3([0,r_min_p(1)],[0,r_min_p(2)],[0,r_min_p(3)],':k'); axis equal
            xlim(2.5*[-1,1]*max(r0_parking_km(:,1)))
            zlim(1.5*[-1,1]*max(r0_parking_km(:,3)))
            hold off
%             view([117,25])
            % view([-123,21])
            view([-90,0])
            label('\itx_L \rm[km]','\ity_L \rm[km]','\itz_L \rm[km]'); title('')
            if thetaf_loop == 99
                % legend([p1,p2,p3,p4,p5,p6],{'主星','初始点','变轨点','目标点','转移轨迹','参考轨迹'},'location','eastoutside')
                legend([p1,p2,p3,p71,p4,p5],{'主星','初始位置','末端位置','最短距离位置','掠飞轨迹','停泊轨迹'},'location','eastoutside')
            end
            % grid on; grid minor; 
            box on
%             title({['{\it\theta}_0=', num2str(theta00), ' rad, {\it\Deltat}=',num2str(theta20f_all(thetaf_loop)), '天']}); 
%             title({['{\it\theta}_0=3\pi/4',' rad, {\it\Deltat}=',num2str(dt_all_day(dt_loop)), '天']}); 
            title({['{\it\theta}_{2,0}^{ 0}=\pi',' rad, {\it\theta}_{2,0}^{ \itf}=\pi rad']}); 
            pause(0.1)
            exportgraphics(gcf,'三维掠飞编队.jpg','Resolution',600)
        end
    end
end

%% 画图
figure(619);

[theta200_p,theta20f_p] = ndgrid(theta200_all,theta20f_all);
surf(theta200_p,theta20f_p,r_min_normal_all);
xlim([min(theta200_all),max(theta200_all)]); ylim([min(theta20f_all),max(theta20f_all)]);
xlabel('{\it\theta}_{2,0}^{ 0} [rad]'); ylabel('{\it\theta}_{2,0}^{ \itf} [rad]'); 
xticks(0:pi/2:2*pi); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks(0:pi/2:2*pi); yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
zlabel('{\ith_z} [km]')
view([42,32])
% plot3(dtheta0_tar, dtheta1_tar, mean(dist_desired),'.r','MarkerSize',20)
% legend('minimum distance','maximum distance','chosen angle difference','Location','north')
% set(gcf,'Color',[255,255,255]/255);
exportgraphics(gcf,'三维掠飞法向距离.jpg','Resolution',600)


figure(621);
surf_dv = surf(theta200_p,theta20f_p,dv_total_m_all);
% set(gca,'ZScale','log')
xlim([min(theta200_all),max(theta200_all)]); ylim([min(theta20f_all),max(theta20f_all)]);
xlabel('{\it\theta}_{2,0}^{ 0} [rad]'); ylabel('{\it\theta}_{2,0}^{ \itf} [rad]'); 
xticks(0:pi/2:2*pi); xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks(0:pi/2:2*pi); yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
zlabel('总脉冲 [m/s]')
view([132,66])
exportgraphics(gcf,'三维掠飞脉冲.jpg','Resolution',600)

