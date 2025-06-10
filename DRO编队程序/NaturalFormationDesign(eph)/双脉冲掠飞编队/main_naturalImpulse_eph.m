clear
clc

format longg
format compact
warning on
dbstop if error

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 13);%坐标轴
set(0,'defaultTextFontSize', 13);%文字
set(0,'defaultLineLineWidth',1.5)

addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

CoreNum = 30; %设定机器CPU核心数量
if isempty(gcp('nocreate')) %如果并行未开启
    parpool(CoreNum);
end
%% 常数与变量
load('generalSolFFT_12.mat')
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);
iDRO = 0;
% load('DROephMultiRev2024.mat')
switch iDRO
case 0
    % DRO1
    load('DROephMultiRev2024.mat')
    x0j2k = x0j2k';
case 1
    % DRO1 2023.6-2027.1 <<1 hr
%     t0UTC = datetime('2023-06-14 07:36:12');
    jd0 = 2460109.81760632; % TBD
    x0j2k = [260861.357166097,167404.871364748,78987.9106263715,...
        -0.797706314601184,0.975263733684134,0.530154969412835];
case 2
    % DRO2 2023.6-2027.6  ~2.5 hr
%     t0UTC = datetime('2023-06-13 23:00:16');
    jd0 = 2460109.45931928;
    x0j2k = [276577.520222802,144514.098050027,67426.0628386977,...
        -0.708023138517267,1.02870635679717,0.558158622414164];
otherwise
    error('no exist DRO')
end

%% CR3BP下的DRO绝对运动
% DRO绝对运动
x0_DRO_CRTBP = para.x0_DRO;
tspan_CRTBP = [0 2*para.T0];
sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),tspan_CRTBP, x0_DRO_CRTBP, opts_ode);
t_sample_CRTBP = linspace(0,para.T0,2000);
xx_MCR_tar_CRTBP = deval(sol_DRO,t_sample_CRTBP)';

isplot = 1; isDisplay = 0;
t_total_eph = 2*para.T0*con.T_norma;

for ii_loop = 1:2
    if ii_loop == 1
        k00 = 5e-5; k0f = -5e-5; 
        theta200 = 0; theta20f = 20/180*pi;
        k2 = 1e-5;
        t01 = 5/con.T_norma_day;
        dt1 = 3/con.T_norma_day;
        tf1 = t01+dt1;
        dt1_parking =6/con.T_norma_day;
        t02 = tf1 + dt1_parking;
    else
        k00 = -5e-5; k0f = 5e-5; 
        theta200 = 20/180*pi; theta20f = pi/2;
        k2 = 1e-5;
        t01 = 14/con.T_norma_day;
        dt1 = 3/con.T_norma_day;
        tf1 = t01+dt1;
        dt1_parking = 10/con.T_norma_day;
        t02 = tf1 + dt1_parking;
    end
    
    xt0f_DRO = deval(sol_DRO,[t01,tf1,t02])';
    
    x01 = generalSol_relMotion(t01,k00,0,k2,0,theta200,para,coe);
    xf1 = generalSol_relMotion(tf1,k0f,0,k2,0,theta20f,para,coe);
    x0_DRO = deval(sol_DRO,t01);
    rv_chaser_km = [x01(1:3)*con.r_norma,xf1(1:3)*con.r_norma; 
    x01(4:6)*con.v_norma,xf1(4:6)*con.v_norma]';
    [x0_DRO_all,r0_REL_all_km,~,r_REL_traj_all_km] = forcedRelMotion(x0_DRO',rv_chaser_km,dt1,'LVLH',0,0);
    dv_all = r0_REL_all_km(:,4:6) - rv_chaser_km(:,4:6);
    dv_all(2,:) = -dv_all(2,:);
    dvnorm_all_km = sqrt(sum(dv_all.^2,2));
    dv_total_m1 = sum(dvnorm_all_km)*1e3;
    
    r_norm_REL_traj_all_km = sqrt(sum(r_REL_traj_all_km(:,1:2).^2,2));
    [r_min_d,r_min_index] = min(r_norm_REL_traj_all_km);
    r_min_p = r_REL_traj_all_km(r_min_index,:);

    t_sample1 = linspace(tf1,t02,200);
    xf1_parking = generalSol_relMotion(t_sample1,k0f,0,k2,0,theta20f,para,coe);
    rf1_parking_km = xf1_parking(1:3,:)'*con.r_norma;
    %% 星历下的主星轨道(DRO绝对运动) 以及 星历下的双脉冲轨道
    aux = []; % 加载星历、设置初始历元
    load('DE430Coeff.mat');%星历表
    aux.C_Mat = DE430Coeff;

    % 星历积分初值（旋转系）
    aux.jd0 = jd0;
    % aux.t0UTC = datetime(2024,1,1,0,0,0);
    aux = initialize(aux); % 初始化

    x0_j2k_tar0 = x0j2k;
    x0_MCR_tar0 = T_ECJ2k2Rot(aux.jd0, x0_j2k_tar0, [0,0,0],aux.C_Mat, 'MCEMR');

    % tspan_eph0 = [0,t_total]; % 10*para.T0*con.T_norma
    tspan_eph0 = [0,t_total_eph]; % 10*para.T0*con.T_norma
    length_t_eph0 = max(ceil(diff(tspan_eph0)/3600)+1,2000);
    t_sample_eph0 = linspace(tspan_eph0(1),tspan_eph0(2),length_t_eph0);
    % 通过主星位置(相位角)确定变轨时间
    theta_tar_eph = mod(atan2(-xt0f_DRO(:,2),xt0f_DRO(:,1)),2*pi);
    % 通过积分event将所有点的相位角等于target_theta的星历DRO上的点标注出来，并得到其时间戳
    MaxStep = para.T0/33*con.T_norma;% 防止错过变轨点
    [~,~,aux_MCR_event1] = Propagate_EphRotFrame(x0_MCR_tar0,tspan_eph0,...
        tspan_eph0,aux,0,@(t, y)DRO_eph_event2(t, y, aux, mod(theta_tar_eph(1)-0.1,2*pi)),MaxStep);
    tspan_eph0 = [aux_MCR_event1.te,t_total_eph];
    [~,~,aux_MCR_event] = Propagate_EphRotFrame(aux_MCR_event1.xe_MCR,tspan_eph0,...
        tspan_eph0,aux,0,@(t, y)DRO_eph_event2(t, y, aux, theta_tar_eph),MaxStep);
    
    % 如果aux_MCR_event存在连续的两个事件,之间的时间差很小，说明是event函数的问题，将第二个事件删掉
    if any(diff(aux_MCR_event.te)<min(para.T0/2*con.T_norma/100,10))
        index_event = find(diff(aux_MCR_event.te)<dt*con.T_norma/100)+1;
        aux_MCR_event.te(index_event) = [];
        aux_MCR_event.ie(index_event) = [];
        aux_MCR_event.xe_j2k(index_event,:) = [];
        aux_MCR_event.ae_j2k(index_event,:) = [];
        aux_MCR_event.xe_MCR(index_event,:) = [];
        aux_MCR_event.ae_MCR(index_event,:) = [];
    end
    
    % 重新积分，得到从变轨起始点至末端点的轨道
    tspan_eph = aux_MCR_event.te([1,end]); % 10*para.T0*con.T_norma
    length_t_eph = max(ceil(diff(tspan_eph)/3600)+1,2000);
    t_sample_eph = linspace(tspan_eph(1),tspan_eph(2),length_t_eph);
    [xx_MCR_tar,aa_MCR_tar] = Propagate_EphRotFrame(aux_MCR_event.xe_MCR(1,:),tspan_eph,...
        t_sample_eph,aux,0);
    % xx_j2k_tar = T_Rot2ECJ2k(aux.jd0+t_sample_eph/86400, xx_MCR_tar, aux.C_Mat, 'MCEMR');
    % xx_ECR_tar = T_ECJ2k2Rot(aux.jd0+t_sample_eph/86400, xx_j2k_tar, [],aux.C_Mat, 'ECEMR');

    % 相位角等于theta_tar_eph(1)的点处，其副星位置设为r1=xx_MCRLVLH_rel_CRTBP_mane_km(1,1:3)'
    % 相位角等于theta_tar_eph(2)的点处，其副星位置设为r2=xx_MCRLVLH_rel_CRTBP_mane_km(2,1:3)'
    r0f_MCRLVLH_rel_eph = [x01(1:3)*con.r_norma,xf1(1:3)*con.r_norma,rf1_parking_km(end,1:3)']';
    % 将第一个变轨点作为新的起始位置
    aux2 = aux;  aux2.t0UTC = []; aux2.jd0 = aux.jd0 + aux_MCR_event.te(1)/86400;
    aux2 = initialize(aux2);
    
    x0_MCR_tar2 = aux_MCR_event.xe_MCR(1,:); % 第一个变轨点处的主星状态
    dt_mane = diff(aux_MCR_event.te);
    
    % 优化
    [xx_MCR_tar_eph_mane,xx_MCRLVLH_rel_eph_mane,dv_all_eph_km,xx_MCR_tar_eph_traj,...
        xx_MCRLVLH_rel_eph_traj,r_MCRLVLH_rel_eph_mane_error,t_sample_traj] = ...
        forcedRelMotion_eph(x0_MCR_tar2,r0f_MCRLVLH_rel_eph,dt_mane(1:2),aux2,'LVLH',0,isDisplay);
    dv_all_eph_m = dv_all_eph_km*1e3;
    dvnorm_all_eph_m = sqrt(sum(dv_all_eph_m.^2,2));

    % 星历的地心惯性系下的LVLH
    xx_j2k_tar_traj = T_Rot2ECJ2k(aux.jd0+t_sample_traj/86400,xx_MCR_tar_eph_traj,aux.C_Mat, 'MCEMR',0);
    aa_j2k_tar_traj = eomj2kMtx(xx_j2k_tar_traj,t_sample_traj,aux);
    xx_MCR_chaser_eph_traj = xx_MCR_tar_eph_traj + ...
        T_TCO2TCR_eph(xx_MCRLVLH_rel_eph_traj,xx_MCR_tar_eph_traj,[],'LVLH');
    xx_j2k_chaser_traj = T_Rot2ECJ2k(aux.jd0+t_sample_traj/86400,xx_MCR_chaser_eph_traj,aux.C_Mat, 'MCEMR',0);
    xx_j2k_rel_traj = xx_j2k_chaser_traj-xx_j2k_tar_traj;
    xx_j2kLVLH_rel_eph_traj = T_TCR2TCO_eph(xx_j2k_rel_traj,xx_j2k_tar_traj,aa_j2k_tar_traj,'LVLH');
    
    t_mane_all = aux_MCR_event.te(1:3) - aux_MCR_event.te(1);
    xx_j2k_tar_mane = T_Rot2ECJ2k(aux.jd0+t_mane_all/86400,xx_MCR_tar_eph_mane,aux.C_Mat, 'MCEMR',0);
    aa_j2k_tar_mane = eomj2kMtx(xx_j2k_tar_mane,t_mane_all,aux);
    xx_MCR_chaser_eph_mane = xx_MCR_tar_eph_mane + ...
        T_TCO2TCR_eph(xx_MCRLVLH_rel_eph_mane,xx_MCR_tar_eph_mane,[],'LVLH');
    xx_j2k_chaser_mane = T_Rot2ECJ2k(aux.jd0+t_mane_all/86400,xx_MCR_chaser_eph_mane,aux.C_Mat, 'MCEMR',0);
    xx_j2k_rel_mane = xx_j2k_chaser_mane-xx_j2k_tar_mane;
    xx_j2kLVLH_rel_eph_mane = T_TCR2TCO_eph(xx_j2k_rel_mane,xx_j2k_tar_mane,aa_j2k_tar_mane,'LVLH');
    %% 画图
    if isplot == 1
        
        f1 = figure(1);
        % xx_MCR_tar = xx_DRO_traj_eph;
        p12 = plot3(xx_MCR_tar(:,1),xx_MCR_tar(:,2),xx_MCR_tar(:,3),'color',[237, 177, 32]/255,'LineWidth',1.5); hold on;
        p13 = plot3(xx_MCR_tar_eph_mane(:,1),xx_MCR_tar_eph_mane(:,2),max(xx_MCR_tar(:,3))+xx_MCR_tar_eph_mane(:,3),'x','color',[54, 100, 208]/255,'LineWidth',1.5);
        p11 = plot3(xx_MCR_tar_CRTBP(:,1)*aux.LU,xx_MCR_tar_CRTBP(:,2)*aux.LU,max(xx_MCR_tar(:,3))+xx_MCR_tar_CRTBP(:,3)*aux.LU,'k','LineWidth',1.0);

        box on; grid on; hold off;
        axis equal;xlabel('\itx_M \rm[km]'); ylabel('\ity_M \rm[km]'); zlabel('\itz_M \rm[km]')
%         L = legend([p11 p12,p13],'CRTBP','星历','变轨位置','Location','northeastoutside');
        L = legend([p11 p12,p13],'CRTBP','Ephemeris','Maneuver Point','Location','northeastoutside');
        % set(L,'box','off')
        xlim([-1.1,1.1]*1e5)
        ylim([-1.3,1.3]*1e5)
        set(gca,'FontSize',15);
        title('DRO')
        % title('DRO (坐标系\itM)')
        % title('DRO (M Frame)')
        set(gcf,'Renderer','painters')
        view(0,90)
        set(gcf,'Color',[255,255,255]/255);
%         export_fig DRO.png -r600
        
        % ---------------------------MCR LVLH坐标系 相对运动------------------------
        % 画图
        figure(2)
        p21 = plot3(0,0,0,'ks');  hold on
        zmax = 0.009;
%         zmax = max(xx_MCRLVLH_rel_eph_traj(:,3));
%         p22 = plot3(xx_MCRLVLH_rel_CRTBP_1_km(:,1),xx_MCRLVLH_rel_CRTBP_1_km(:,2),zmax+xx_MCRLVLH_rel_CRTBP_1_km(:,3),'k','LineWidth',1);
%         plot3(xx_MCRLVLH_rel_CRTBP_2_km(:,1),xx_MCRLVLH_rel_CRTBP_2_km(:,2),zmax+xx_MCRLVLH_rel_CRTBP_2_km(:,3),'k','LineWidth',1);
        p24 = plot3(xx_MCRLVLH_rel_eph_traj(:,1),xx_MCRLVLH_rel_eph_traj(:,2),xx_MCRLVLH_rel_eph_traj(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5);
        p25 = plot3(xx_MCRLVLH_rel_eph_mane(1:2,1),xx_MCRLVLH_rel_eph_mane(1:2,2),xx_MCRLVLH_rel_eph_mane(1:2,3),'x','color',[54, 100, 208]/255,'LineWidth',1.5); 

%         set(gcf, 'Units', 'centimeters', 'Position',[10 10 7 5]);
%         L = legend([p21 p22,p24,p25],{'主星','CRTBP','星历','变轨点'});
%         L.Location = 'northeastoutside'; L.AutoUpdate = 'off';
%         L = legend([p21 p22,p24,p25],{'Origin','CRTBP','Ephemeris','maneuver point'}); 
%         L.Orientation = 'horizontal';
%         grid on; 
        axis equal; box on; 
%         hold off
%         xlim(scale*[-0.801,0.801]); ylim(scale*[-1,1]); 
%         xlim(scale*[-0.901,0.901]); ylim(scale*[-1.301,1.301]); 
%         title('相对运动')
%         title(['\Delta{\itt}=',num2str(dt/para.T0),'{\itT}_0',])
%         view([0,90])
%         view(-55,30)
        set(gca,'FontSize',15);
%         set(gcf,'Units','centimeters','Position',[35,15,11,9]); % 2D formation
        xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); zlabel('\itz_L \rm[km]');
%         xlabel('\itx_L \rm[km]','Units','centimeters','Position',[6.4,-0.35,0]); % 3D formation
%         set(gcf,'Units','centimeters','Position',[35,15,11,7.3]); % 3D formation
        
%         set(gcf,'Color',[255,255,255]/255);
%         export_fig 2DForma.png -r600

        % ---------------------------j2kLVLH坐标系 相对运动-----------------------------
        % 画图
%         figure(3)
%         p21 = plot3(0,0,0,'ks');  hold on
% %         zmax = 0.009*scale;
%         p22 = plot3(rr_ECILVLH_rel_CRTBP_1(:,1),rr_ECILVLH_rel_CRTBP_1(:,2),zmax+rr_ECILVLH_rel_CRTBP_1(:,3),'k','LineWidth',1);
%         plot3(rr_ECILVLH_rel_CRTBP_2(:,1),rr_ECILVLH_rel_CRTBP_2(:,2),zmax+rr_ECILVLH_rel_CRTBP_2(:,3),'k','LineWidth',1);
%         p24 = plot3(xx_j2kLVLH_rel_eph_traj(:,1),xx_j2kLVLH_rel_eph_traj(:,2),xx_j2kLVLH_rel_eph_traj(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5);
%         p25 = plot3(xx_j2kLVLH_rel_eph_mane(:,1),xx_j2kLVLH_rel_eph_mane(:,2),xx_j2kLVLH_rel_eph_mane(:,3),'x','color',[54, 100, 208]/255,'LineWidth',1.5); 
% 
%         set(gcf, 'Units', 'centimeters', 'Position',[23 6 15 10]);
%         L = legend([p21 p22,p24,p25],{'主航天器','CRTBP','星历','变轨点'});
%         L.Location = 'northeastoutside'; L.AutoUpdate = 'off';
%         grid on; axis equal; box on; hold off
% %         xlim(scale*[-0.801,0.801]); ylim(scale*[-1,1]); 
% %         xlim(scale*[-0.501,0.501]); ylim(scale*[-1.201,0.101]); 
% %         title('相对运动')
% %         title(['\Delta{\itt}=',num2str(dt/para.T0),'{\itT}_0',])
% %         view([0,90])
%         view(-57,20)
%         set(gca,'FontSize',15);
% %         set(gcf,'Units','centimeters','Position',[35,15,11,9]); % 2D formation
%         xlabel('\itx_L_E \rm[km]'); ylabel('\ity_L_E \rm[km]'); zlabel('\itz_L_E \rm[km]');
% 
%         set(gcf,'Color',[255,255,255]/255);
% %         export_fig 2DFormaJ2kVLH.png -r600
    end
end
% clear DE430Coeff aux aux2
% save(FileName)


