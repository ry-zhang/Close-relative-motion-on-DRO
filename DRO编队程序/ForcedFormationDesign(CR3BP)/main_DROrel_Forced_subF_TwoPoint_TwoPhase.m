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

sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 3*para.T0], para.x0_DRO, opts);

%% 选择patch pointd 的方式
k1 = 0; theta1_0 = 0.625*pi; 
k2 = 10e-6;

for ii_loop = 1:2
    if ii_loop == 1
        k00 = 5e-5; k0f = -5e-5; 
        theta200 = 0; theta20f = 20/180*pi;
        t01 = 5/con.T_norma_day;
        dt1 = 3/con.T_norma_day;
        tf1 = t01+dt1;
        dt1_parking =6/con.T_norma_day;
        t02 = tf1 + dt1_parking;
    else
        k00 = -5e-5; k0f = 5e-5; 
        theta200 = 20/180*pi; theta20f = pi/2;
        t01 = 14/con.T_norma_day;
        dt1 = 3/con.T_norma_day;
        tf1 = t01+dt1;
        dt1_parking = 10/con.T_norma_day;
        t02 = tf1 + dt1_parking;
    end
    isplot = 1;



    x01 = generalSol_relMotion(t01,k00,0,k2,0,theta200,para,coe);
    xf1 = generalSol_relMotion(tf1,k0f,0,k2,0,theta20f,para,coe);
    x0_DRO = deval(sol_DRO,t01);
    rv_chaser_km = [x01(1:3)*con.r_norma,xf1(1:3)*con.r_norma; 
    x01(4:6)*con.v_norma,xf1(4:6)*con.v_norma]';
    [x0_DRO_all,r0_REL_all_km,~,r_REL_traj_all_km] = forcedRelMotion(x0_DRO',rv_chaser_km,dt1,Flag_coordinate,0,0);
    dv_all = r0_REL_all_km(:,4:6) - rv_chaser_km(:,4:6);
    dv_all(2,:) = -dv_all(2,:);
    dvnorm_all_km = sqrt(sum(dv_all.^2,2));
    dv_total_m1 = sum(dvnorm_all_km)*1e3;
    
    r_norm_REL_traj_all_km = sqrt(sum(r_REL_traj_all_km(:,1:2).^2,2));
    [r_min_d,r_min_index] = min(r_norm_REL_traj_all_km);
    r_min_p = r_REL_traj_all_km(r_min_index,:);

    t_sample1 = linspace(tf1,t02,200);
    xf1_parking = generalSol_relMotion(t_sample1,k0f,0,k2,0,theta20f,para,coe);


    if isplot == 1
    %     t_sample2 = linspace(tf2,tf2+dt2_parking,200);
    %     xf2_parking = generalSol_relMotion(t_sample2,k00,0,k2,0,theta20f,para,coe);
        rf1_parking_km = xf1_parking(1:3,:)'*con.r_norma;
    %     rf2_parking_km = xf2_parking(1:3,:)'*con.r_norma;
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
        p5 = plot3(rf1_parking_km(:,1),rf1_parking_km(:,2),rf1_parking_km(:,3),'Color',colorTable(1,:));
    %     p6 = plot3(rf2_parking_km(:,1),rf2_parking_km(:,2),rf2_parking_km(:,3),'Color',colorTable(1,:)); axis equal
        p71 = plot3(r_min_p(1),r_min_p(2),r_min_p(3),'ok'); axis equal
        p72 = plot3([0,r_min_p(1)],[0,r_min_p(2)],[0,r_min_p(3)],':k'); 
        axis equal
%         xlim(2.5*[-1,1]*max(rf1_parking_km(:,1)))
%         zlim(1.5*[-1,1]*max(rf1_parking_km(:,3)))
    %     hold off
    %             view([117,25])
        % view([-123,21])
        view([-90,0])
        label('\itx_L \rm[km]','\ity_L \rm[km]','\itz_L \rm[km]'); title('')

            % legend([p1,p2,p3,p4,p5,p6],{'主星','初始点','变轨点','目标点','转移轨迹','参考轨迹'},'location','eastoutside')
%             legend([p1,p2,p3,p71,p4,p5],{'主星','初始位置','末端位置','最短距离位置','掠飞轨迹','停泊轨迹'},'location','eastoutside')

        % grid on; grid minor; 
        box on
    %             title({['{\it\theta}_0=', num2str(theta00), ' rad, {\it\Deltat}=',num2str(theta20f_all(thetaf_loop)), '天']}); 
    %             title({['{\it\theta}_0=3\pi/4',' rad, {\it\Deltat}=',num2str(dt_all_day(dt_loop)), '天']}); 
%         title({['{\it\theta}_{2,0}^{ 0}=\pi',' rad, {\it\theta}_{2,0}^{ \itf}=\pi rad']}); 
        pause(0.1)

    end
end
hold off
exportgraphics(gcf,'两段掠飞编队.jpg','Resolution',600)