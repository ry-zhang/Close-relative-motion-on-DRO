% 给定ECILVLH下的相位alpha，输出DRO在MCR下的相位theta
% 定义alpha为B星在ECILVLH中的相位，从ECILVLH的x轴正向逆时针旋转的角度，范围为0~360
% 定义theta为A星在MCR中的相位，从MCR的x轴正向逆时针旋转的角度，范围为0~360

function theta_tar = alpha2theta(alpha_tar)

%% MCR下的DRO前后跟飞编队，也即CRTBP下的MCR DRO周期相对运动
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% --------------------加载CRTBP下的轨道数据-------------------------
load('DRO1to2_FloquetEigen')
% load('FloquetEig_all.mat'); ii_DRO = 71;
x0_DRO = x0_DRO_M_3d;
x0_REL_pos = EigenVector.p3'/con.r_norma/EigenVector.p3(2); % 只考虑一个方向的周期相对运动

% ----------------------------在MCR与MCRLVLH下积分----------------------------
T0 = para.T0;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 T0], [x0_DRO, x0_REL_pos], opts_ode);
t_sample_CRTBP = linspace(0,T0,2000);
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

%% 拟合ECI下的相对运动参考周期轨道
% ----------------------------拟合自变量与因变量--------------------------
% theta为主星在MCR坐标系中的几何相位(逆时针)
theta_all = mod(atan2(-xx_MCR_target_CRTBP(:,2),xx_MCR_target_CRTBP(:,1)),2*pi);
% alpha为副星在ECILVLH坐标系中的几何相位(逆时针)
alpha_all = mod(atan2(-xx_ECILVLH_rel_CRTBP(:,2),xx_ECILVLH_rel_CRTBP(:,1)),2*pi);

% 滑动theta与alpha，方便插值
theta_all_shift = mod(theta_all - pi,2*pi);
theta_all_shift(end) = 2*pi;
alpha_all_shift = mod(alpha_all - 3/2*pi,2*pi);
alpha_all_shift(end) = 2*pi;

%% alpha2theta
r2d = 180/pi;

alpha_tar_shift = mod(alpha_tar/r2d - 3/2*pi,2*pi);
theta_tar_shift = interp1(alpha_all_shift,theta_all_shift,alpha_tar_shift);
theta_tar = mod(theta_tar_shift + pi,2*pi)*r2d;

% theta_tar_shift = mod(theta_tar/r2d - pi,2*pi);
% alpha_tar_shift = interp1(theta_all_shift,alpha_all_shift,theta_tar_shift);
% alpha_tar = mod(alpha_tar_shift + 3/2*pi,2*pi)*r2d;



