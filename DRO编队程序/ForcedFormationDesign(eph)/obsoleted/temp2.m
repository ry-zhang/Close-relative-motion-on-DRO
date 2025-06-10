clear
addpath('../subF_eom(CR3BP)')
addpath('../subF_eom(eph)')

format longg
format compact
warning off

%% 常数与变量
load('FloquetEig')
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);
opts_fsolve = optimoptions('fsolve','Display','iter',...
    'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'Algorithm','levenberg-marquardt');

%% CR3BP下的DRO周期相对运动
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = -2/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';

% 周期轨道
tf = para.T0;
sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
length_t = 2000;
t_sample = linspace(0,tf,length_t);
PeriodicTraj = deval(sol,t_sample)';
rDRO_norma = PeriodicTraj(:,1:3);
rRelMotion_norma = PeriodicTraj(:,7:9);

% 星历下的月球相对于地球位置
aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC  = [2030 1 1 0 0 0]; % 初始历元
aux = initialize(aux); % 初始化
% 计算
t_sample_day = t_sample*con.T_norma_day;
tt_jd = aux.jd0 + t_sample_day;
rEarth_all = zeros(length_t,3);
rMoon_all = zeros(length_t,3);
for ii = 1:length(tt_jd)
    t_jd_ii = tt_jd(ii);
    pEarthRV  = JPL_Eph_DE430_PosVel(t_jd_ii,  3, aux.C_Mat); % 地球
    rEarth_all(ii,:) = pEarthRV(1:3);
    pMoonRV  = JPL_Eph_DE430_PosVel(t_jd_ii,  10, aux.C_Mat); % 月球
    rMoon_all(ii,:) = pMoonRV(1:3);
end
% rMEbarycener = (rEarth_all*aux.xmu(3)+rMoon_all*aux.xmu(10))/(aux.xmu(3)+aux.xmu(10));
rMoon2Ear = rMoon_all - rEarth_all;
rMoon2Ear_norm = sqrt(sum(rMoon2Ear.^2,2));

%% 采用星历下月球距地球的距离转换DRO周期相对运动
rRelMotion_km_norma = rRelMotion_norma.*con.r_norma;
rRelMotion_km = rRelMotion_norma.*rMoon2Ear_norm;

%% 画图
figure(1)
plot(rRelMotion_km_norma(:,1), rRelMotion_km_norma(:,2),'k','LineWidth',2); hold on;

plot(rRelMotion_km(:,1), rRelMotion_km(:,2),'color',[0 0.4470 0.7410]); hold on;
plot(rRelMotion_km(1,1), rRelMotion_km(1,2),'g^');
plot(rRelMotion_km(end,1), rRelMotion_km(end,2),'rv');
box on; grid on; grid minor; hold off; 
axis equal; xlabel('x[km]'); ylabel('y[km]')
set(gca,'FontSize',15,'fontname','times new roman'); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title('RelMotion(LVLH)')

