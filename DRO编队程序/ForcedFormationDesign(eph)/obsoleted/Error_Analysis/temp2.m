clear
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

dt = 10*3600;

% 初始化
aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC  = [2030 1 1 0 0 0]; % 初始历元
aux = initialize(aux); % 初始化
t0UTC = aux.t0UTC;

x0_MCR = [72687.2175909459 0 0 0 -0.540957166578942 0];
% x0_MCR = [7000,0,0,0,7,0];
x0_j2k = T_Rot2ECJ2k(aux.jd0,x0_MCR,aux.C_Mat,'MCEMR');

% 从x0_j2k递推一段,得到t2时刻的状态x0_MCR_target2
t2 = dt;
tspan_sec = [0,t2];
x0_j2k2 = Propagate_EphJ2kFrame(x0_j2k,tspan_sec,t2,aux);

% 从x0_j2k递推更长时间，得到t2:t3之间的轨道
t3 = t2+20*dt;
t_sample = linspace(t2,t3,2000);
xx_j2k_tall = Propagate_EphJ2kFrame(x0_j2k,[0,t3],t_sample,aux);

% 从x0_j2k2递推至t3，同样的，得到t2:t3之间的轨道
aux.jd0 = aux.jd0 + t2/86400;
tspan_sec2 = [0,t3-t2];
t_sample2 = t_sample-t2;
xx_j2k_loop = Propagate_EphJ2kFrame(x0_j2k2,tspan_sec2,t_sample2,aux);

% 求两次积分之间的误差
error_j2k = xx_j2k_tall - xx_j2k_loop;
error_j2k_norm = sqrt(sum(error_j2k(:,1:3).^2,2));
plot(t_sample2/3600,error_j2k_norm)
xlabel('t [hr]'); ylabel('error [km]')