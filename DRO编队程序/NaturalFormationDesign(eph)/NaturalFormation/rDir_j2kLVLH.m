function phi_nagetive = rDir_j2kLVLH(t0,x0_MCR_target,auxSFF)
% 给定t0、x0_MCRLVLH_target，可求出r0_MCRLVLH_rel_CRTBP,即可得到x0_j2kLVLH_rel的方向，即直线方向
% 在其方向上逆时针偏转13.5°，即可得到phi_nagetive，分离脉冲沿这个方向，后续变轨总脉冲会位于最小变轨脉冲附近
% 根据计算结果，可确定初始时刻优化方向
% 因为求phi_nagetive的时候，是按照离轨时刻的主星状态计算的，加上最大1天的变轨时间，会逆时针偏转[0,27]°左右
% 注意j2kLVLH下的相对运动周期与DRO周期是相同的，也即为13.6天左右

aux = []; 
aux.jd0 = auxSFF.jd0;
aux = initialize(aux);
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;

%% CR3BP下的DRO周期相对运动
load('FloquetEig12')
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = -1/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';

% 周期轨道
tf = para.T0;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
xx_MCR_target_CRTBP = sol_deval_CRTBP(:,1:6);
rr_MCRLVLH_rel_CRTBP = sol_deval_CRTBP(:,7:9);

%% 拟合参考周期轨道
rr_MCRLVLH_rel_CRTBP_km = rr_MCRLVLH_rel_CRTBP.*con.r_norma;
theta_all = mod(atan2(-xx_MCR_target_CRTBP(:,2),xx_MCR_target_CRTBP(:,1)),2*pi);

% Set up fittype and options.
ft = fittype( 'fourier8' );
opts_fit = fitoptions( 'Method', 'NonlinearLeastSquares','Display','Off' );

% Fit model to data.
fitresult_x = fit( theta_all, rr_MCRLVLH_rel_CRTBP_km(:,1), ft, opts_fit );
fitresult_y = fit( theta_all, rr_MCRLVLH_rel_CRTBP_km(:,2), ft, opts_fit );

%% 计算对应的x0_j2k_target参考轨道点 
x0_j2k_target = T_Rot2ECJ2k(aux.jd0+t0/86400,x0_MCR_target,aux.C_Mat, 'MCEMR',0);
theta_MCR_target = mod(atan2(-x0_MCR_target(:,2),x0_MCR_target(:,1)),2*pi);
fit_x = feval(fitresult_x,theta_MCR_target);
fit_y = feval(fitresult_y,theta_MCR_target);
x0_MCRLVLH_rel = [fit_x,fit_y,zeros(1,4)];
x0_MCR_chaser = T_TCO2TCR_eph(x0_MCRLVLH_rel,x0_MCR_target,[],'LVLH') + x0_MCR_target;
x0_j2k_chaser = T_Rot2ECJ2k(aux.jd0+t0/86400,x0_MCR_chaser,aux.C_Mat, 'MCEMR',0);
x0_j2k_rel = x0_j2k_chaser-x0_j2k_target;
x0_j2kLVLH_rel = T_TCR2TCO_eph(x0_j2k_rel,x0_j2k_target,[],'LVLH');
phi_nagetive = atan2d(x0_j2kLVLH_rel(2),x0_j2kLVLH_rel(1)) - 13.5;

