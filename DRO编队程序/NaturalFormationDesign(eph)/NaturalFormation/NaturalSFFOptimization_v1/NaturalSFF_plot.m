clear
% close all
addpath('../../../subF_eom(CR3BP)')
addpath('../../../subF_eom(eph)')

format longg; format compact
warning off

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

load('FloquetEig12')
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

%% 画图定义值
load('NaturalSSF_DRO1_100day_scale90_SEMephemeris.mat')

num = 1;
scale_0 = 1; % 初值放缩
% numDay = 365;
% tspan_sec = [0,86400*numDay];
% t_sample = linspace(tspan_sec(1),tspan_sec(2),1500);

%% CR3BP下的DRO周期相对运动
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = -1/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';

% 周期轨道
tf = para.T0;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
rvDRO_MCR_ref = sol_deval_CRTBP(:,1:6);

rRelMot_ref_L = sol_deval_CRTBP(:,7:9);
rRelMot_ref = rRelMot_ref_L(:,[2,3,1]);
rRelMot_ref(:,[2,3]) = -rRelMot_ref(:,[2,3]);

vRelMot_ref_L = sol_deval_CRTBP(:,10:12);
vRelMot_ref = vRelMot_ref_L(:,[2,3,1]);
vRelMot_ref(:,[2,3]) = -vRelMot_ref(:,[2,3]);

rRelMot_km_L = rRelMot_ref_L.*con.r_norma;

%% 计算优化结果
aux = []; 
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;

aux.t0UTC = t0UTC;
if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,1); % 初始化
else
    aux = initialize(aux); % 初始化
end

% 将DRO初值向后递推epoch_SSF秒
epoch_SSF = epoch_SSF_all(num);
tspan_sec_re0 = tspan_sec+epoch_SSF;
t_sample_re0 = t_sample+epoch_SSF;

if epoch_SSF == 0
    xx_j2k_target0 = x0_j2k_target0;
else
    if strcmp(DynamicModel,'HighFidelity')
        [~,xx_j2k_target0] = ode113(@eom_SPICE_stm , [0,epoch_SSF], x0_j2k_target0 , opts_ode, aux);
    else
        [~,xx_j2k_target0] = ode113(@eqom_geoMEMEJ2k, [0,epoch_SSF], x0_j2k_target0 , opts_ode, aux);
    end
end
x0_j2k_target = xx_j2k_target0(end,:);

% 初始化新的aux
% aux = aux; 
% aux.t0UTC = [];
% aux.jd0 = aux.jd0 + epoch_SSF/86400;
% aux = initialize(aux); % 初始化

% 主星星历积分
x0_MCR_target = T_ECJ2k2Rot(aux.jd0 + epoch_SSF/86400, x0_j2k_target, [0,0,0],aux.C_Mat, 'MCEMR');
[xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec_re0,t_sample_re0,aux);

x0_TCO_rel_optimal = scale_0*x0_TCO_optimize_all(num,:);
x0_TCR_chaser = [x0_TCO_rel_optimal(1),0,x0_TCO_rel_optimal(2),x0_TCO_rel_optimal(3),0,x0_TCO_rel_optimal(4)];
x0_MCR_chaser = T_TCO2TCR_eph(x0_TCR_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec_re0,t_sample_re0,aux);
% 计算相对运动
rv_MCR_rel = xx_MCR_chaser-xx_MCR_target;
rv_TCO_rel = T_TCR2TCO_eph(rv_MCR_rel,xx_MCR_target,a_MCR_target,'LVLH');

%% LVLH坐标系下的相对运动
%---------月心旋转系下的LVLH-----------
figure(1)
plot3(0,0,0,'ks','MarkerSize',5); hold on;
rRelMot_km_L_s = scale_0*scale_ref*rRelMot_km_L;
plot3(rRelMot_km_L_s(:,1), rRelMot_km_L_s(:,2), rRelMot_km_L_s(:,3),'k','LineWidth',1.5); 
plot3(rv_TCO_rel(:,1), rv_TCO_rel(:,2), rv_TCO_rel(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5);
plot3(rv_TCO_rel(1,1), rv_TCO_rel(1,2), rv_TCO_rel(1,3),'g^');
plot3(rv_TCO_rel(end,1), rv_TCO_rel(end,2), rv_TCO_rel(end,3),'rv');
legend('主星','CRTBP','星历','初始位置','终端位置','Location','northeastoutside');
box on; grid on; grid minor; hold off; 
axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]')
ylim(scale_ref*[-1.3,0.2]); xlim(scale_ref*[-1,1]);
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title('月心旋转系下的LVLH坐标系')
view(0,90)

%---------地心惯性系下的LVLH-----------
figure(2)
% CRTBP的地心惯性系下的LVLH
rvDRO_ECR_CRTBP = [sol_deval_CRTBP(:,1:3)'-[-1,0,0]'; sol_deval_CRTBP(:,4:6)'];
rvDRO_ECI_CRTBP = synodic2inertial(rvDRO_ECR_CRTBP,t_sample_CRTBP)';
rRelMot_ECR_linear = T_TCO2TCR_CR3BP(rRelMot_km_L_s,rvDRO_MCR_ref,'LVLH',con.mu)'; 
rRelMot_ECI = synodic2inertial(rRelMot_ECR_linear,t_sample_CRTBP)';
rv_TCO_ECI_f = T_TCI2TCO_E_CR3BP(rRelMot_ECI,rvDRO_ECI_CRTBP,t_sample_CRTBP,'LVLH',con.mu)';
% 星历的地心惯性系下的LVLH
if strcmp(DynamicModel,'HighFidelity')
    sol_j2k_target = ode113(@eom_SPICE_stm , tspan_sec_re0, x0_j2k_target , opts_ode, aux);
else
    sol_j2k_target = ode113(@eqom_geoMEMEJ2k, tspan_sec_re0, x0_j2k_target , opts_ode, aux);
end

[xx_j2k_target,a_j2k_target] = deval(sol_j2k_target,t_sample_re0);
xx_j2k_target = xx_j2k_target'; a_j2k_target = a_j2k_target(1:3,:)';
% xx_j2k_target = T_Rot2ECJ2k(aux.jd0+(t_sample_re0)/86400,xx_MCR_target,aux.C_Mat, 'MCEMR',0);
xx_j2k_chaser = T_Rot2ECJ2k(aux.jd0+t_sample_re0/86400,xx_MCR_chaser,aux.C_Mat, 'MCEMR',0);
rv_RelMot_j2k = xx_j2k_chaser-xx_j2k_target;
rv_RelMot_j2k_LVLH = T_TCR2TCO_eph(rv_RelMot_j2k,xx_j2k_target,a_j2k_target,'LVLH');
rv_RelMot_j2k_LVLH = rv_RelMot_j2k_LVLH';

% 画图
plot3(0,0,0,'ks','MarkerSize',5); hold on
plot3(rv_TCO_ECI_f(1,:),rv_TCO_ECI_f(2,:),rv_TCO_ECI_f(3,:),'k','LineWidth',1.5);
plot3(rv_RelMot_j2k_LVLH(1,:),rv_RelMot_j2k_LVLH(2,:),rv_RelMot_j2k_LVLH(3,:),'LineWidth',1.5,'Color',[237, 177, 32]/255);
plot3(rv_RelMot_j2k_LVLH(1,1),rv_RelMot_j2k_LVLH(2,1),rv_RelMot_j2k_LVLH(3,1),'g^');
plot3(rv_RelMot_j2k_LVLH(1,end),rv_RelMot_j2k_LVLH(2,end),rv_RelMot_j2k_LVLH(3,end),'rv');
hold off; 
box on; grid on; grid minor; 
axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]')
legend('主星','CRTBP','星历','初始位置','终端位置','Location','northeastoutside');
title('地心J2000下的LVLH坐标系')
set(gca,'FontSize',15);

%% DRO 绝对运动
%---------月心旋转系下的DRO-----------
f3 = figure(3);
p2 = plot3(xx_MCR_target(:,1),xx_MCR_target(:,2),xx_MCR_target(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5); hold on;
p3 = plot3(xx_MCR_target(1,1),xx_MCR_target(1,2),xx_MCR_target(1,3),'g^');
p4 = plot3(xx_MCR_target(end,1),xx_MCR_target(end,2),xx_MCR_target(end,3),'rv');
p1 = plot3(rvDRO_MCR_ref(:,1)*con.r_norma,rvDRO_MCR_ref(:,2)*con.r_norma,rvDRO_MCR_ref(:,3)*con.r_norma,'k','LineWidth',1.5);
box on; grid on; grid minor; hold off;
axis equal;xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]'); zlabel('\itz \rm[km]')
% L = legend([p1 p2,p3,p4],'CRTBP','星历','初始位置','终端位置','Location','northeast');
L = legend([p1 p2,p3,p4],'CRTBP','星历','初始位置','终端位置','Location','northeastoutside');
% set(L,'box','off')
xlim([-1.2,1.2]*1e5)
set(gca,'FontSize',15);
title('DRO (MCR)')
% title('DRO (M Frame)')
set(gcf,'Renderer','painters')
view(0,90)

%---------月心旋转系下的DRO-----------
xx_j2k_target = T_Rot2ECJ2k(aux.jd0+t_sample_re0/86400,xx_MCR_target,aux.C_Mat, 'MCEMR',0);
t_sample_temp = linspace(0,2*para.T0,2000);
dphase = -para.T0/4.8;
[~,rDRO_MCR_ref] = ode113(@(t,x)eom_abs3b(t,x,con.mu),t_sample_temp, x0_DRO_M_3d, opts_ode); % 注意这个旋转系与MCR不一样

% 星历下DRO与CRTBP下DRO在MCR坐标系下是相近的，但是在转化到ECI的时候，
% 星历下考虑的地月距离波动较大，而CRTBP下考虑的地月距离是恒定的，
% 因此二者在ECI下会差别过大，除非将CRTBP转化的时候采用实时的地月距离。
% 由于LVLH坐标系下是已经优化过的，因此在LVLH坐标系下的相对轨迹反而差别较小

rDRO_ECI_ref = synodic2inertial((rDRO_MCR_ref(:,1:3)'-[-1,0,0]'),dphase+t_sample_temp)'*con.r_norma;
f4 = figure(4);
set(f4,'name','星历DRO ECI')
% figure('color',[1 1 1],'name','星历DRO MCR')
p2 = plot3(xx_j2k_target(:,1),xx_j2k_target(:,2),xx_j2k_target(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5); hold on;
p3 = plot3(xx_j2k_target(1,1),xx_j2k_target(1,2),xx_j2k_target(1,3),'g^');
p4 = plot3(xx_j2k_target(end,1),xx_j2k_target(end,2),xx_j2k_target(end,3),'rv');
% p1 = plot3(rDRO_ECI_ref(:,1),rDRO_ECI_ref(:,2),rDRO_ECI_ref(:,3),'k','LineWidth',1.5);
box on; grid on; grid minor; hold off;
axis equal;xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]'); zlabel('\itz \rm[km]')
L = legend([p2,p3,p4],'CRTBP','星历','初始位置','终端位置','Location','northeastoutside');
% set(L,'box','off')
% xlim([-1.2,1.2]*1e5)
set(gca,'FontSize',15);
title('DRO (J2000)')
set(gcf,'Renderer','painters')
view(0,90)
