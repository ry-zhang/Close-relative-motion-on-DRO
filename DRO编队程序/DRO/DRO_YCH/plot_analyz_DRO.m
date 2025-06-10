clear
close all
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
load('DRO_all3.mat')
T_Moon = 27.2844292118811;
J_all = J_period_all(:,1);
T_all = J_period_all(:,2)*T_Moon;
r_norma = 384399;

%% 能量与周期
figure(1)
plot(J_all,T_all,'LineWidth',1.5); hold on
xlabel('雅各比积分常数','FontSize',13);
ylabel('\itT \rm[day]')
set(gca,'FontSize',13);
grid on
xlim([min(J_all),max(J_all)]);
ylim([min(T_all),max(T_all)]);
plot([min(J_all),max(J_all)],[T_Moon/2,T_Moon/2],'--','LineWidth',1.5,'Color',[217, 83, 25]/255)
plot([min(J_all),max(J_all)],[T_Moon/3,T_Moon/3],'--','LineWidth',1.5,'Color',[217, 83, 25]/255)
hold off

%% 能量与近月点/远月点
global mu
mu = 0.01215;
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);
r_peri_all = zeros(size(J_all));
r_apo_all = zeros(size(J_all));
for ii_loop = 1:1:length(delta_T_all)
    DT = T_DRO_all(ii_loop)*2*pi;
    t_sample = linspace(0,DT,300);
    state_ini = state_ini_all(ii_loop,:);
    sol = ode113(@(t,x)pcr3bp(t,x),[0 DT], state_ini, opts_ode); 
    sol_sample = deval(sol,t_sample);
    r_DRO = sol_sample(1:2,:)'-[1-mu,0]; % 月心旋转系
%     r_DRO = sol_sample(1:2,:)'; % 共同质心旋转系
    r_DRO_norm = sqrt(sum(r_DRO.^2,2));
    r_peri_all(ii_loop) = min(r_DRO_norm);
    r_apo_all(ii_loop) = max(r_DRO_norm);
end
figure(2)
plot(J_all,r_peri_all*r_norma,'LineWidth',1.5); hold on
plot(J_all,r_apo_all*r_norma,'LineWidth',1.5); 
xlabel('雅各比积分常数','FontSize',13);
ylabel('距离 \rm[km]')
set(gca,'FontSize',13);
grid on
xlim([min(J_all),max(J_all)]);
legend('近月点','远月点')
