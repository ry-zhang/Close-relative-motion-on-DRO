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
% load('FloquetEig_all.mat')
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);


%% CR3BP下的DRO周期相对运动
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = -1/con.r_norma/EigenVector.p3(2)*EigenVector.p3';

tf = J_period_all(4,2)*2*pi*10;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
xx_MCR_target_CRTBP = sol_deval_CRTBP(:,1:6);
rr_MCRLVLH_rel_CRTBP = sol_deval_CRTBP(:,7:9);

% CRTBP的地心惯性系下的DRO
t_sample_temp = linspace(0,tf,2000);
dphase = 0;
[~,rv_MCR_target] = ode113(@(t,x)eom_abs3b(t,x,con.mu),t_sample_temp, x0_DRO_loop, opts_ode); % 注意这个旋转系与MCR不一样
rv_ECI_target = synodic2inertial(([rv_MCR_target(:,1:3)'-[-1,0,0]';rv_MCR_target(:,4:6)']),t_sample_temp)';
% r0_ECI_chaser = synodic2inertial((rr_MCR_target(2,1:3)'-[-1,0,0]'),t_sample_temp(1))'*con.r_norma;
rv0_MCR_chaser = inertial2synodic((rv_ECI_target(2,:)'),t_sample_temp(1))'+[-1,0,0,0,0,0];
[~,rv_MCR_chaser] = ode113(@(t,x)eom_abs3b(t,x,con.mu),t_sample_temp, rv0_MCR_chaser, opts_ode); % 注意这个旋转系与MCR不一样
rv_MCR_rel = rv_MCR_chaser - rv_MCR_target;
rv_MCRLVLH_rel = T_TCR2TCO_CR3BP(rv_MCR_rel,rv_MCR_target,'LVLH',con.mu);

figure(1)
plot3(rv_MCRLVLH_rel(:,1),rv_MCRLVLH_rel(:,2),rv_MCRLVLH_rel(:,3));
axis equal
title('MCRLVLH')
view(0,90)

figure(2)
plot3(rv_MCR_target(:,1),rv_MCR_target(:,2),rv_MCR_target(:,3),'k'); hold on
plot3(rv_MCR_chaser(:,1),rv_MCR_chaser(:,2),rv_MCR_chaser(:,3)); hold off
axis equal
title('MCR')
view(0,90)
