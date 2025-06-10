%% linear relative motion error analysis near DRO orbit
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../../subF_eom(CR3BP)')
% 加载2:1共振DRO的线性化相对运动
load('FloquetEig12.mat')

%% definite constants and variables
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2
mu = 0.01215; % 
r_norma = 3.84399*10^5; % km
% T_M = 27.321661; % day
T_norma = sqrt(r_norma^3/(mu_E+mu_M)); % s
T_norma_day = T_norma/3600/24;
v_norma = r_norma/T_norma; % km/s
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
opts_event = odeset('RelTol',1e-13,'AbsTol',1e-20,'Events',@eventsPeriod);

%% 数值积分
% initial value of rel motion in LVLH frame
% Sol_linear为线性化相对运动，vec1-2分别是周期解、发散解，vec3-4是平面拟周期、vec5-6是法向拟周期
r0_REL = Sol_linear.vec3(1:3)' * 1e-5; 
v0_REL = Sol_linear.vec3(4:6)' * 1e-5; 
x0_REL = [r0_REL,v0_REL];

% 线性化相对运动积分
T0 = para.T0;
dt = 10*T0; % integration time
% dt = 3*T0; % integration time
% integration
sol = ode113(@(t,x)eom_rel3b(t,x,mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts_event);

length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*T_norma_day;
sol_sample = deval(sol,t_sample);
abs_motion_M = sol_sample(1:6,:);
rel_motion_L_linear = sol_sample(7:12,:);


% 将相同的相对运动初值转换至旋转系中，得到副星准确的绝对运动
x1_DRO_M_3d = x0_DRO_M_3d + T_TCR2TCO_CR3BP([r0_REL,v0_REL],sol_sample(1:6,1)','LVLH',con.mu); 
sol2 = ode113(@(t,x)eom_abs3b(t,x,mu),[0 dt], x1_DRO_M_3d, opts_event);
sol2_sample = deval(sol2,t_sample);
rel_motion_M_Abs = sol2_sample(1:6,:)-sol_sample(1:6,:);
% 绝对运动作差转移至LVLH坐标系
rel_motion_L_Abs = T_TCO2TCR_CR3BP(rel_motion_M_Abs',abs_motion_M','LVLH',con.mu)';

% 线性化相对运动的误差（归一化单位）
error = rel_motion_L_Abs - rel_motion_L_linear;

