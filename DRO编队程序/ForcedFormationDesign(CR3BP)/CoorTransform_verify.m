%% 验证坐标变换的正确性
clear
close all

%% 常数
para.mu = 0.01215; % 20200531
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2
para.r_norma = 3.84399*10^5; % km
para.T_norma = sqrt(para.r_norma^3/(mu_E+mu_M)); % s
para.T_norma_day = para.T_norma/3600/24;
para.v_norma = para.r_norma/para.T_norma; % km/s

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 绝对运动与相对运动初值
x0_target = [-0.191872994233853,0,0,0,0.531458007318788,0]; % DRO
T0 = 3.39476129234148;
% x0_target = [0.162164232531584,0,0.0791678086999886,0,-0.192755627370783,0]; % Halo
% T0 = 2.77600473787838;
% 积分轨道,计算单值矩阵中心流形（拟周期项）
sol_MOT = ode113(@(t,x)eomM_abs3b(t,x,para.mu),[0 1*T0], [x0_target,reshape(eye(6),1,36)],opts_ode);
M0T = reshape(sol_MOT.y(7:42,end),6,6); % 单值矩阵
[Mveig,Meig] = eig(M0T);
x0_rel_M = real(Mveig(:,3))'*-2/para.r_norma;
x0_rel_L = T_TCR2TCO_CR3BP(x0_rel_M,sol_MOT.y(1:6,end)','LVLH',para.mu);

dt = T0;
sol = ode113(@(t,x)eom_rel3b(t,x,para.mu),[0 dt], [x0_target,x0_rel_L], opts_ode);
t_sample = linspace(0,dt,20000);
sol_deval = deval(sol,t_sample)';

%% 坐标转换
x_DRO_all = sol_deval(:,1:6);

rv_rel_L = [sol_deval(:,7:9),sol_deval(:,10:12)];
rv_rel_M = T_TCO2TCR_CR3BP(rv_rel_L,x_DRO_all,'LVLH',para.mu);
% rv_rel_L = T_TCR2TCO_CR3BP(rv_rel_M,x_DRO_all,'LVLH',para.mu);
rv_rel_V = T_TCR2TCO_CR3BP(rv_rel_M,x_DRO_all,'VVLH',para.mu);
rv_rel_C = T_TCR2TCO_CR3BP(rv_rel_M,x_DRO_all,'VNC',para.mu);
% rv_rel_L2 = T_TCO2TCR(rv_rel_M,x_DRO_all,'LVLH',para.mu);
% norm(rv_rel_L2-rv_rel_L)

rv_TCO_c = T_TCO2TCR_CR3BP(rv_rel_L,x_DRO_all,'LVLH',para.mu);
% rv_TCO_c2 = T_TCR2TCO(rv_rel_L,x_DRO_all,'LVLH',para.mu);
% norm(rv_TCO_c2-rv_TCO_c)


%% 画图
set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultTextFontSize', 14)
% set(0, 'DefaultLineLineWidth', 1)
set(0, 'DefaultAxesBox', 'on')
% set(0, 'DefaultAxesDataAspectRatio', [1,1,1]) % axis equal
set(0, 'DefaultAxesXGrid', 'on'); set(0, 'DefaultAxesYGrid', 'on'); set(0, 'DefaultAxesZGrid', 'on')

figure(1)
r_DRO_km = sol_deval(:,1:3)*para.r_norma;
plot3(r_DRO_km(:,1),r_DRO_km(:,2),r_DRO_km(:,3),'Color',[0, 114, 189]/255,'LineWidth',1.5); 
axis equal; box on; view(0,90)
xlabel('x_M [km]'); ylabel('y_M [km]'); zlabel('z_M [km]')

figure(2)
subplot(2,2,1)
plot3(0,0,0,'ks');  hold on
plot3(rv_rel_M(:,1)*para.r_norma,rv_rel_M(:,2)*para.r_norma,rv_rel_M(:,3)*para.r_norma,'Color',[0, 114, 189]/255,'LineWidth',1.5); 
axis equal; box on; view(0,90)
xlabel('x_M [km]'); ylabel('y_M [km]'); zlabel('z_M [km]')
subplot(2,2,2)
plot3(0,0,0,'ks');  hold on
plot3(rv_rel_L(:,1)*para.r_norma,rv_rel_L(:,2)*para.r_norma,rv_rel_L(:,3)*para.r_norma,'Color',[0, 114, 189]/255,'LineWidth',1.5); 
axis equal; box on; view(0,90)
xlabel('x_L [km]'); ylabel('y_L [km]'); zlabel('z_L [km]')
subplot(2,2,3)
plot3(0,0,0,'ks');  hold on
plot3(rv_rel_V(:,1)*para.r_norma,rv_rel_V(:,2)*para.r_norma,rv_rel_V(:,3)*para.r_norma,'Color',[0, 114, 189]/255,'LineWidth',1.5); 
axis equal; box on; view(0,0)
xlabel('x_V [km]'); ylabel('y_V [km]'); zlabel('z_V [km]')
subplot(2,2,4)
plot3(0,0,0,'ks');  hold on
plot3(rv_rel_C(:,1)*para.r_norma,rv_rel_C(:,2)*para.r_norma,rv_rel_C(:,3)*para.r_norma,'Color',[0, 114, 189]/255,'LineWidth',1.5); 
axis equal; box on; view(0,0)
xlabel('x_C [km]'); ylabel('y_C [km]'); zlabel('z_C [km]')

figure(3)
var_tested = rv_rel_M;
R = diff(var_tested(2:end,1:3))./var_tested(2:end-1,4:6)/diff(t_sample(1:2));
plot(t_sample(3:end)'*para.T_norma_day,R,'.','LineWidth',1.5); 
xlabel('t [day]')
xlim([0,max(t_sample)*para.T_norma_day]); ylim([0.9,1.1])
legend('R_x','R_y','R_z')