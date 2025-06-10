% 程序中航天器状态变量命名规则
% -------------------------------------------------------------------------
% 命名格式: 变量名_坐标系_航天器_后缀
% 变量名:   x0-初始状态，xf-末端状态，xx-整个轨道状态，对于a/r/v，采用同样命名规则
% 坐标系:   MCR-月心旋转系，j2k-地心J2000坐标系，MCRLVLH-MCR下的LVLH相对系，j2kLVLH-j2k下的LVLH相对系
% 航天器:   target-DRO目标航天器，chaser-绕飞航天器，rel-两者相对运动
% 后缀:     可用于标识额外信息，可有多段。例如描述不同阶段：Last\Next\transfer,于描述不同模型：CRTBP/eph
% 例:       xx_MCR_chaser_Next \ x0_MCRLVLH_rel_Last \ xf_j2kLVLH_rel_Last

dbstop if error

clear
% close all
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg; format compact
warning off

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

load('FloquetEig12')
% load('FloquetEig_all.mat'); ii_DRO = 71;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

%% 画图定义值
load('PhDNaturalSSF_DRO0_SEMephemeris_posi.mat')
naturalSFF = auxSFF.naturalSFF;
ii_segment = 1;
% numDay = 530;
% tspan_sec = naturalSFF(ii_segment).t0+[0,86400*numDay];

%% CR3BP下的DRO周期相对运动
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = auxSFF.flag_refDirec/con.r_norma/EigenVector.p3(2)*EigenVector.p3';

tf = J_period_all(4,2)*2*pi;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
xx_MCR_target_CRTBP = sol_deval_CRTBP(:,1:6);
rr_MCRLVLH_rel_CRTBP = sol_deval_CRTBP(:,7:9);

%% 计算优化结果
aux = []; 
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.jd0 = auxSFF.jd0;
if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,1); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux = initialize(aux); % 初始化
else
    error('Wrong DynamicModel')
end

% 将DRO初值向后递推epoch_SSF秒
scale_ref = naturalSFF(ii_segment).scale_ref;
% 优化时间段
epoch_SSF = naturalSFF(ii_segment).t0;
tspan_sec = [epoch_SSF,naturalSFF(ii_segment).tf];
t_sample = linspace(tspan_sec(1),tspan_sec(2),max(500,ceil(diff(tspan_sec)/3600/4)+1));

x0_MCR_target = auxSFF.naturalSFF(ii_segment).x0_MCR_target;
[xx_MCR_target,aa_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);

x0_MCRLVLH_rel = naturalSFF(ii_segment).x0_MCRLVLH_rel;
x0_MCR_chaser = T_TCO2TCR_eph(x0_MCRLVLH_rel,x0_MCR_target,aa_MCR_target(1,:),'LVLH')+x0_MCR_target;
% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);
% 计算相对运动
xx_MCR_rel = xx_MCR_chaser-xx_MCR_target;
xx_MCRLVLH_rel = T_TCR2TCO_eph(xx_MCR_rel,xx_MCR_target,aa_MCR_target,'LVLH');

%% 计算轨道数据
% 月心旋转系下的LVLH
rr_MCRLVLH_rel_CRTBP_km = scale_ref*rr_MCRLVLH_rel_CRTBP.*con.r_norma;

% CRTBP的地心惯性系下的LVLH
xx_ECR_target_CRTBP = [sol_deval_CRTBP(:,1:3)'-[-1,0,0]'; sol_deval_CRTBP(:,4:6)'];
xx_ECI_target_CRTBP = synodic2inertial(xx_ECR_target_CRTBP,t_sample_CRTBP)';
rr_MCR_rel_CRTBP = T_TCO2TCR_CR3BP(rr_MCRLVLH_rel_CRTBP_km,xx_MCR_target_CRTBP,'LVLH',con.mu)'; 
rr_ECI_rel_CRTBP = synodic2inertial(rr_MCR_rel_CRTBP,t_sample_CRTBP)';
rr_ECILVLH_rel_CRTBP = T_TCI2TCO_E_CR3BP(rr_ECI_rel_CRTBP,xx_ECI_target_CRTBP,t_sample_CRTBP,'LVLH',con.mu);

% 星历的地心惯性系下的LVLH
xx_j2k_target = T_Rot2ECJ2k(aux.jd0+t_sample/86400,xx_MCR_target,aux.C_Mat, 'MCEMR',0);
aa_j2k_target = eomj2kMtx(xx_j2k_target,t_sample,aux);
xx_j2k_chaser = T_Rot2ECJ2k(aux.jd0+t_sample/86400,xx_MCR_chaser,aux.C_Mat, 'MCEMR',0);
xx_j2k_rel = xx_j2k_chaser-xx_j2k_target;
xx_j2kLVLH_rel = T_TCR2TCO_eph(xx_j2k_rel,xx_j2k_target,aa_j2k_target,'LVLH');

%% LVLH坐标系下的相对运动
%---------月心旋转系下的LVLH-----------
figure(1)
plot3(0,0,0,'ks','MarkerSize',5); hold on;
plot3(rr_MCRLVLH_rel_CRTBP_km(:,1), rr_MCRLVLH_rel_CRTBP_km(:,2), max(xx_MCRLVLH_rel(:,3))+rr_MCRLVLH_rel_CRTBP_km(:,3),'k','LineWidth',1.5); 
p2 = plot3(xx_MCRLVLH_rel(:,1), xx_MCRLVLH_rel(:,2), xx_MCRLVLH_rel(:,3),'Color',[237, 177, 32]/255,'LineWidth',1);
% p2.Color(4) = 0.4;
plot3(xx_MCRLVLH_rel(1,1), xx_MCRLVLH_rel(1,2), max(xx_MCRLVLH_rel(:,3))+xx_MCRLVLH_rel(1,3),'g^');
plot3(xx_MCRLVLH_rel(end,1), xx_MCRLVLH_rel(end,2), max(xx_MCRLVLH_rel(:,3))+xx_MCRLVLH_rel(end,3),'rv');
box on; grid on; grid minor; 
hold off; 
axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); zlabel('\itz_L \rm[km]')
yLimit = scale_ref*((auxSFF.flag_refDirec == -1)*[-1.3,0.2] + (auxSFF.flag_refDirec == 1)*[-0.2,1.3]);
ylim(yLimit); 
% xlim(scale_ref*[-1,1]);
% title('星历下编队(MCR LVLH)')
xlim(scale_ref*[-0.5,0.5]);
switch ii_segment
    case 1
        title('伴飞编队1')
    case 2
        title('伴飞编队2')
    case 3
        title('伴飞编队3')
        legend('主星','CRTBP','星历','初始位置','终端位置','Location','northeastoutside');
end
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
view(0,90)
% set(gcf,'Color',[255,255,255]/255);
% export_fig MCRLVLH.png -r600

%---------地心惯性系下的LVLH-----------
figure(2)
plot3(0,0,0,'ks','MarkerSize',5); hold on;
plot3(rr_ECILVLH_rel_CRTBP(:,1), rr_ECILVLH_rel_CRTBP(:,2), max(xx_j2kLVLH_rel(:,3))+rr_ECILVLH_rel_CRTBP(:,3),'k','LineWidth',1.5); 
p2 = plot3(xx_j2kLVLH_rel(:,1), xx_j2kLVLH_rel(:,2), xx_j2kLVLH_rel(:,3),'Color',[237, 177, 32]/255,'LineWidth',1);
% p2.Color(4) = 0.4;
plot3(xx_j2kLVLH_rel(1,1), xx_j2kLVLH_rel(1,2), max(xx_j2kLVLH_rel(:,3))+xx_j2kLVLH_rel(1,3),'g^');
plot3(xx_j2kLVLH_rel(end,1), xx_j2kLVLH_rel(end,2), max(xx_j2kLVLH_rel(:,3))+xx_j2kLVLH_rel(end,3),'rv');
box on; grid on; grid minor; 
hold off; 
axis equal; xlabel('\itx_{LE} \rm[km]'); ylabel('\ity_{LE} \rm[km]'); zlabel('\itz_{LE} \rm[km]')
yLimit = scale_ref*[-1.3,1.3];
ylim(yLimit); 
% xlim(scale_ref*[-1,1]);
% title('星历下编队(MCR LVLH)')
xlim(scale_ref*[-0.83,0.83]);
switch ii_segment
    case 1
        title('伴飞编队1')
    case 2
        title('伴飞编队2')
    case 3
        title('伴飞编队3')
        legend('主星','CRTBP','星历','初始位置','终端位置','Location','northeastoutside');
end
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
view(0,90)
set(gcf,'Color',[255,255,255]/255);
export_fig ECILVLH.png -r600


%% DRO 绝对运动
%---------月心旋转系下的DRO-----------
figure(3);
p2 = plot3(xx_MCR_target(:,1),xx_MCR_target(:,2),xx_MCR_target(:,3),'Color',[237, 177, 32]/255,'LineWidth',1); hold on;
p2.Color(4) = 0.6;
p3 = plot3(xx_MCR_target(1,1),xx_MCR_target(1,2),xx_MCR_target(1,3),'g^');
p4 = plot3(xx_MCR_target(end,1),xx_MCR_target(end,2),xx_MCR_target(end,3),'rv');
p1 = plot3(xx_MCR_target_CRTBP(:,1)*con.r_norma,xx_MCR_target_CRTBP(:,2)*con.r_norma,max(xx_MCR_target(:,3))+xx_MCR_target_CRTBP(:,3)*con.r_norma,'k','LineWidth',1.5);
box on; grid on; grid minor; 
hold off;
axis equal;xlabel('\itx_M \rm[km]'); ylabel('\ity_M \rm[km]'); zlabel('\itz_M \rm[km]')
L = legend([p1 p2,p3,p4],'CRTBP','星历','初始位置','终端位置','Location','northeastoutside');
% set(L,'box','off')
xlim([-0.9,0.9]*1e5)
set(gca,'FontSize',13);
title('DRO (MCR)')
% title('DRO (M Frame)')
% set(gcf,'Renderer','painters')
view(0,90)
% view(-42,28)
% set(gcf,'Color',[255,255,255]/255);
% export_fig DRO.png -r600

