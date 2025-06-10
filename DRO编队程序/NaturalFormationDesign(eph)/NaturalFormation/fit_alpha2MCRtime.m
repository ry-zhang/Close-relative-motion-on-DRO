% 给定ECILVLH下的相位，输出DRO在MCR下的相位(时间)[0,T]

function t_tar = fit_alpha2MCRtime(alpha_tar)

%% MCR下的DRO前后跟飞编队，也即CRTBP下的MCR DRO周期相对运动
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% --------------------加载CRTBP下的轨道数据-------------------------
load('FloquetEig12')
% load('FloquetEig_all.mat'); ii_DRO = 71;
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop_pos = EigenVector.p3'/con.r_norma/EigenVector.p3(2);

% ----------------------------在MCR与MCRLVLH下积分----------------------------
tf = J_period_all(4,2)*2*pi;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop_pos], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
xx_MCR_target_CRTBP = sol_deval_CRTBP(:,1:6);
xx_MCRLVLH_rel_CRTBP = sol_deval_CRTBP(:,7:12);

% ----------------------------主星轨道转移至ECI----------------------------
xx_ECR_target_CRTBP = xx_MCR_target_CRTBP;
xx_ECR_target_CRTBP(:,1:3) = xx_MCR_target_CRTBP(:,1:3) - [-1,0,0];
xx_ECI_target_CRTBP = synodic2inertial(xx_ECR_target_CRTBP',t_sample_CRTBP)';

% ----------------------------相对运动轨道转移至ECILVLH----------------------------
xx_MCR_rel_CRTBP = T_TCO2TCR_CR3BP(xx_MCRLVLH_rel_CRTBP,xx_MCR_target_CRTBP,'LVLH',con.mu);
xx_ECI_rel_CRTBP = synodic2inertial(xx_MCR_rel_CRTBP',t_sample_CRTBP)';
xx_ECILVLH_rel_CRTBP = T_TCI2TCO_E_CR3BP(xx_ECI_rel_CRTBP,xx_ECI_target_CRTBP,t_sample_CRTBP,'LVLH',con.mu);
rr_ECILVLH_rel_CRTBP = xx_ECILVLH_rel_CRTBP(:,1:3);
% vv_ECILVLH_rel_CRTBP = vv_ECILVLH_rel_CRTBP(:,4:6);

%% 拟合ECI下的相对运动参考周期轨道
% ----------------------------拟合自变量与因变量--------------------------
% theta为主星在MCR坐标系中的几何相位(逆时针)
% theta_all = mod(atan2(-xx_MCR_target_CRTBP(:,2),xx_MCR_target_CRTBP(:,1)),2*pi);
% alpha为副星在ECILVLH坐标系中的几何相位(顺时针)
alpha_all = mod(atan2(-rr_ECILVLH_rel_CRTBP(:,2),rr_ECILVLH_rel_CRTBP(:,1)),2*pi);
alpha_all_shift = mod(alpha_all - 3/2*pi,2*pi);
alpha_all_shift(end) = 2*pi;
% -----------------------------------拟合---------------------------------
% ft = fittype( 'fourier8' );
% opts_fit = fitoptions( 'Method', 'NonlinearLeastSquares','Display','Off' );
% [fit_alpha2t, gof_alpha2t] = fit( alpha_all_shift,t_sample_CRTBP'*con.T_norma_day , ft, opts_fit );


%% 计算初始状态在ECILVLH下的几何相位，以及ECILVLH下的目标几何相位的时间
r2d = 180/pi;
alpha_tar_shift = mod(alpha_tar/r2d - 3/2*pi,2*pi);
% t_tar = feval(fit_alpha2t,alpha_tar/r2d);
t_tar = interp1(alpha_all_shift,t_sample_CRTBP'*con.T_norma_day,alpha_tar_shift);
% alpha_f_deg = mod(alpha_f_shift + 3/2*pi,2*pi)*r2d;




