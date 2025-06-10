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

CoreNum =6; %设定机器CPU核心数量
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
tspan_CRTBP = [0 para.T0];
sol_tar_CRTBP = ode113(@(t,x)eom_abs3b(t,x,con.mu),tspan_CRTBP, x0_DRO_CRTBP, opts_ode);
t_sample_CRTBP = linspace(0,para.T0,2000);
xx_MCR_tar_CRTBP = deval(sol_tar_CRTBP,t_sample_CRTBP)';

isplot = 1; isDisplay = 0;
dt_all = linspace(0.001,0.800,800)*para.T0;

% isplot = 1; isDisplay = 1;
% dt_all = [0.001,0.01,0.1,0.2]*para.T0;
% dt_all = [0.3,0.4,0.5,0.6]*para.T0;
% dt_all = [0.3,0.4,0.55,0.58]*para.T0;

t_total_eph = 11.7*para.T0*con.T_norma;

scale_all = [1,10,100,1000]; % km
for jj_index = 3
% for jj_index = 1
    scale = scale_all(jj_index);
    
    dv_CRTBP_all = zeros(length(dt_all),1);
    dv_eph_all = zeros(length(dt_all),24);
    r_MCRLVLH_rel_eph_mane_error_all = dv_eph_all;

% for ii_index = 1:length(dt_all)
% for ii_index = [300,310,320]
% for ii_index = [313,553]
% for ii_index = [100,200,400,700]
for ii_index = 200
    tic
    %% CRTBP中的双脉冲轨道
    % 计算自然轨道
    theta1_0 = 0.625*pi; % 1.625*pi 
    theta2_0 = 0.266*pi; % 1.266*pi 
%     k0 = -0.197; k1 = 1; k2 = 0; Dimflag = '2D';
    k0 = -0.197; k1 = 1; k2 = 0.5; Dimflag = '3D';
    FileName = [Dimflag,'Formation_DR0',num2str(iDRO),'_Scale',num2str(scale)];
    
%     dt = 0.1*para.T0;
    dt = dt_all(ii_index);
    
    t_sample_CRTBP_1 = linspace(dt/2,para.T0-dt/2,200);
    xx_MCRLVLH_rel_CRTBP_1 = generalSol_relMotion(t_sample_CRTBP_1,k0,k1,k2,theta1_0,theta2_0,para,coe)';
    xx_MCR_tar_CRTBP_1 = deval(sol_tar_CRTBP,t_sample_CRTBP_1)';

    % 将y轴最大值放缩为1km
    x0_MCRLVLH_rel_halfT = generalSol_relMotion(para.T0/2,k0,k1,k2,theta1_0,theta2_0,para,coe)';
    ratio = scale/abs(x0_MCRLVLH_rel_halfT(2));
    xx_MCRLVLH_rel_CRTBP_1 = ratio*xx_MCRLVLH_rel_CRTBP_1;
    xx_MCRLVLH_rel_CRTBP_1_km = [xx_MCRLVLH_rel_CRTBP_1(:,1:3),...
        xx_MCRLVLH_rel_CRTBP_1(:,4:6)*con.v_norma/con.r_norma];
    xx_MCR_rel_CRTBP_1 = T_TCO2TCR_CR3BP(xx_MCRLVLH_rel_CRTBP_1,xx_MCR_tar_CRTBP_1,'LVLH',con.mu);
    xx_MCR_rel_CRTBP_1_km = [xx_MCR_rel_CRTBP_1(:,1:3),...
        xx_MCR_rel_CRTBP_1(:,4:6)*con.v_norma/con.r_norma];

    % 计算变轨脉冲
    x0_MCR_tar = deval(sol_tar_CRTBP,para.T0-dt/2)'; % inverse MCR
    [xx_MCR_tar_CRTBP_mane,xx_MCRLVLH_rel_CRTBP_mane_km] = ...
        forcedRelMotion(x0_MCR_tar,xx_MCRLVLH_rel_CRTBP_1_km([end,1],1:3),dt,'LVLH',0,isDisplay);

    % 积分变轨轨道
    t_sample_CRTBP_2 = linspace(para.T0-dt/2,para.T0+dt/2,200);
    x0_MCRLVLH_rel_CRTBP_2 = [xx_MCRLVLH_rel_CRTBP_mane_km(1,1:3),...
        xx_MCRLVLH_rel_CRTBP_mane_km(1,4:6)*con.r_norma/con.v_norma];
    [~,xx_CRTBP_2] = ode113(@(t,x)eom_rel3b(t,x,con.mu),t_sample_CRTBP_2, [x0_MCR_tar,x0_MCRLVLH_rel_CRTBP_2], opts_ode);
    xx_MCR_tar_CRTBP_2 = xx_CRTBP_2(:,1:6);
    xx_MCRLVLH_rel_CRTBP_2 = xx_CRTBP_2(:,7:12);
    xx_MCRLVLH_rel_CRTBP_2_km = [xx_MCRLVLH_rel_CRTBP_2(:,1:3),...
        xx_MCRLVLH_rel_CRTBP_2(:,4:6)*con.v_norma/con.r_norma];
    xx_MCR_rel_CRTBP_2 = T_TCO2TCR_CR3BP(xx_MCRLVLH_rel_CRTBP_2,xx_MCR_tar_CRTBP_2,'LVLH',con.mu);
    xx_MCR_rel_CRTBP_2_km = [xx_MCR_rel_CRTBP_2(:,1:3),...
        xx_MCR_rel_CRTBP_2(:,4:6)*con.v_norma/con.r_norma];

    % 计算两次变轨脉冲大小
    dv_all_CRTBP_m = (xx_MCRLVLH_rel_CRTBP_mane_km(1:2,4:6)-xx_MCRLVLH_rel_CRTBP_1_km([end,1],4:6))*10^3;
    dvnorm_all_CRTBP_m = sqrt(sum(dv_all_CRTBP_m.^2,2));
    dv_CRTBP_all(ii_index) = mean(dvnorm_all_CRTBP_m);
    
    % 将变轨末端点的速度换为自然轨道段的速度，以为星历模型提供更好的初值
    xx_MCRLVLH_rel_CRTBP_mane_km(2,:) = xx_MCRLVLH_rel_CRTBP_1_km(1,:);

    % 将CRTBP轨道数据转换至地心LVLH坐标系
    xx_ECR_tar_CRTBP_1 = [xx_MCR_tar_CRTBP_1(:,1:3)-[-1,0,0], xx_MCR_tar_CRTBP_1(:,4:6)];
    xx_ECI_tar_CRTBP_1 = synodic2inertial(xx_ECR_tar_CRTBP_1',t_sample_CRTBP_1)';
    rr_ECI_rel_CRTBP_1 = synodic2inertial(xx_MCR_rel_CRTBP_1',t_sample_CRTBP_1)';
    rr_ECILVLH_rel_CRTBP_1 = T_TCI2TCO_E_CR3BP(rr_ECI_rel_CRTBP_1,xx_ECI_tar_CRTBP_1,t_sample_CRTBP_1,'LVLH',con.mu);
    xx_ECR_tar_CRTBP_2 = [xx_MCR_tar_CRTBP_2(:,1:3)-[-1,0,0], xx_MCR_tar_CRTBP_2(:,4:6)];
    xx_ECI_tar_CRTBP_2 = synodic2inertial(xx_ECR_tar_CRTBP_2',t_sample_CRTBP_2)';
    rr_ECI_rel_CRTBP_2 = synodic2inertial(xx_MCR_rel_CRTBP_2',t_sample_CRTBP_2)';
    rr_ECILVLH_rel_CRTBP_2 = T_TCI2TCO_E_CR3BP(rr_ECI_rel_CRTBP_2,xx_ECI_tar_CRTBP_2,t_sample_CRTBP_2,'LVLH',con.mu);

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
    theta_tar_eph = mod(atan2(-xx_MCR_tar_CRTBP_mane(:,2),xx_MCR_tar_CRTBP_mane(:,1)),2*pi);
    % 通过积分event将所有点的相位角等于target_theta的星历DRO上的点标注出来，并得到其时间戳
    MaxStep = min(dt,para.T0-dt)*con.T_norma*3/4;% 防止错过变轨点
    [~,~,aux_MCR_event] = Propagate_EphRotFrame(x0_MCR_tar0,tspan_eph0,...
        t_sample_eph0,aux,0,@(t, y)DRO_eph_event(t, y, aux, theta_tar_eph),MaxStep);
    
    % 如果aux_MCR_event存在连续的两个事件,之间的时间差很小，说明是event函数的问题，将第二个事件删掉
    if any(diff(aux_MCR_event.te)<min(dt*con.T_norma/100,10))
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
    r0f_MCRLVLH_rel_eph = ones(size(aux_MCR_event.ie,2),1).*(aux_MCR_event.ie'==1)*xx_MCRLVLH_rel_CRTBP_mane_km(1,:)+...
        ones(size(aux_MCR_event.ie,2),1).*(aux_MCR_event.ie'==2)*xx_MCRLVLH_rel_CRTBP_mane_km(2,:);
    % 将第一个变轨点作为新的起始位置
    aux2 = aux;  aux2.t0UTC = []; aux2.jd0 = aux.jd0 + aux_MCR_event.te(1)/86400;
    aux2 = initialize(aux2);
    
    x0_MCR_tar2 = aux_MCR_event.xe_MCR(1,:); % 第一个变轨点处的主星状态
    dt_mane = diff(aux_MCR_event.te);
    
    % 优化
    [xx_MCR_tar_eph_mane,xx_MCRLVLH_rel_eph_mane,dv_all_eph_km,xx_MCR_tar_eph_traj,...
        xx_MCRLVLH_rel_eph_traj,r_MCRLVLH_rel_eph_mane_error,t_sample_traj] = ...
        forcedRelMotion_eph(x0_MCR_tar2,r0f_MCRLVLH_rel_eph,dt_mane,aux2,'LVLH',0,isDisplay);
    dv_all_eph_m = dv_all_eph_km*1e3;
    dvnorm_all_eph_m = sqrt(sum(dv_all_eph_m.^2,2));
    dv_eph_all(ii_index,1:length(dvnorm_all_eph_m)) = dvnorm_all_eph_m';
    r_MCRLVLH_rel_eph_mane_error_all(ii_index,1:length(r_MCRLVLH_rel_eph_mane_error)) = r_MCRLVLH_rel_eph_mane_error';
    
    % 星历的地心惯性系下的LVLH
    xx_j2k_tar_traj = T_Rot2ECJ2k(aux.jd0+t_sample_traj/86400,xx_MCR_tar_eph_traj,aux.C_Mat, 'MCEMR',0);
    aa_j2k_tar_traj = eomj2kMtx(xx_j2k_tar_traj,t_sample_traj,aux);
    xx_MCR_chaser_eph_traj = xx_MCR_tar_eph_traj + ...
        T_TCO2TCR_eph(xx_MCRLVLH_rel_eph_traj,xx_MCR_tar_eph_traj,[],'LVLH');
    xx_j2k_chaser_traj = T_Rot2ECJ2k(aux.jd0+t_sample_traj/86400,xx_MCR_chaser_eph_traj,aux.C_Mat, 'MCEMR',0);
    xx_j2k_rel_traj = xx_j2k_chaser_traj-xx_j2k_tar_traj;
    xx_j2kLVLH_rel_eph_traj = T_TCR2TCO_eph(xx_j2k_rel_traj,xx_j2k_tar_traj,aa_j2k_tar_traj,'LVLH');
    
    t_mane_all = aux_MCR_event.te - aux_MCR_event.te(1);
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
        zmax = 0.009*scale;
%         zmax = max(xx_MCRLVLH_rel_eph_traj(:,3));
        p22 = plot3(xx_MCRLVLH_rel_CRTBP_1_km(:,1),xx_MCRLVLH_rel_CRTBP_1_km(:,2),zmax+xx_MCRLVLH_rel_CRTBP_1_km(:,3),'k','LineWidth',1);
        plot3(xx_MCRLVLH_rel_CRTBP_2_km(:,1),xx_MCRLVLH_rel_CRTBP_2_km(:,2),zmax+xx_MCRLVLH_rel_CRTBP_2_km(:,3),'k','LineWidth',1);
        p24 = plot3(xx_MCRLVLH_rel_eph_traj(:,1),xx_MCRLVLH_rel_eph_traj(:,2),xx_MCRLVLH_rel_eph_traj(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5);
        p25 = plot3(xx_MCRLVLH_rel_eph_mane(:,1),xx_MCRLVLH_rel_eph_mane(:,2),xx_MCRLVLH_rel_eph_mane(:,3),'x','color',[54, 100, 208]/255,'LineWidth',1.5); 

%         set(gcf, 'Units', 'centimeters', 'Position',[10 10 7 5]);
%         L = legend([p21 p22,p24,p25],{'主星','CRTBP','星历','变轨点'});
%         L.Location = 'northeastoutside'; L.AutoUpdate = 'off';
%         L = legend([p21 p22,p24,p25],{'Origin','CRTBP','Ephemeris','maneuver point'}); 
%         L.Orientation = 'horizontal';
        grid on; axis equal; box on; hold off
%         xlim(scale*[-0.801,0.801]); ylim(scale*[-1,1]); 
        xlim(scale*[-0.901,0.901]); ylim(scale*[-1.301,1.301]); 
%         title('相对运动')
        title(['\Delta{\itt}=',num2str(dt/para.T0),'{\itT}_0',])
        view([0,90])
%         view(-55,30)
        set(gca,'FontSize',15);
        set(gcf,'Units','centimeters','Position',[35,15,11,9]); % 2D formation
        xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); zlabel('\itz_L \rm[km]');
%         xlabel('\itx_L \rm[km]','Units','centimeters','Position',[6.4,-0.35,0]); % 3D formation
%         set(gcf,'Units','centimeters','Position',[35,15,11,7.3]); % 3D formation
        
        set(gcf,'Color',[255,255,255]/255);
        export_fig 2DForma.png -r600

        % ---------------------------j2kLVLH坐标系 相对运动-----------------------------
        % 画图
        figure(3)
        p21 = plot3(0,0,0,'ks');  hold on
        zmax = 0.009*scale;
        p22 = plot3(rr_ECILVLH_rel_CRTBP_1(:,1),rr_ECILVLH_rel_CRTBP_1(:,2),zmax+rr_ECILVLH_rel_CRTBP_1(:,3),'k','LineWidth',1);
        plot3(rr_ECILVLH_rel_CRTBP_2(:,1),rr_ECILVLH_rel_CRTBP_2(:,2),zmax+rr_ECILVLH_rel_CRTBP_2(:,3),'k','LineWidth',1);
        p24 = plot3(xx_j2kLVLH_rel_eph_traj(:,1),xx_j2kLVLH_rel_eph_traj(:,2),xx_j2kLVLH_rel_eph_traj(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5);
        p25 = plot3(xx_j2kLVLH_rel_eph_mane(:,1),xx_j2kLVLH_rel_eph_mane(:,2),xx_j2kLVLH_rel_eph_mane(:,3),'x','color',[54, 100, 208]/255,'LineWidth',1.5); 

        set(gcf, 'Units', 'centimeters', 'Position',[23 6 15 10]);
        L = legend([p21 p22,p24,p25],{'主航天器','CRTBP','星历','变轨点'});
        L.Location = 'northeastoutside'; L.AutoUpdate = 'off';
        grid on; axis equal; box on; hold off
%         xlim(scale*[-0.801,0.801]); ylim(scale*[-1,1]); 
        xlim(scale*[-0.501,0.501]); ylim(scale*[-1.201,0.101]); 
%         title('相对运动')
%         title(['\Delta{\itt}=',num2str(dt/para.T0),'{\itT}_0',])
%         view([0,90])
        view(-57,20)
        set(gca,'FontSize',15);
%         set(gcf,'Units','centimeters','Position',[35,15,11,9]); % 2D formation
        xlabel('\itx_L_E \rm[km]'); ylabel('\ity_L_E \rm[km]'); zlabel('\itz_L_E \rm[km]');

        set(gcf,'Color',[255,255,255]/255);
        export_fig 2DFormaJ2kVLH.png -r600
    end
    t_end = toc;
    fprintf('第%d次优化,耗时%.1f秒，共%d次优化\n',ii_index,t_end,length(dt_all))
end
% clear DE430Coeff aux aux2
% save(FileName)
end

