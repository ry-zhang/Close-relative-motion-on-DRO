clear
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg
format compact
warning off

%% 定轨误差分析
aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC  = [2023 1 1 0 0 0]; % 初始历元
aux = initialize(aux); % 初始化
t0UTC = aux.t0UTC;

FrameFlag = 'VVLH';

load DRO1_LiAISON.dat 
load DRO1_real.dat 
load DRO2_LiAISON.dat 
load DRO2_real.dat
load LEO_real.dat

size_Data = size(DRO2_real,1);
data_sample = [1:10:size_Data];
t_sample = data_sample*60; % 秒，采样点
size_Sample = length(t_sample);
tt_jd = aux.jd0 + t_sample/86400;

rv_DRO1_real = DRO1_real(data_sample,7:12);
rv_DRO2_real = DRO2_real(data_sample,7:12);
rv_DRO1_LiAISON = DRO1_LiAISON(data_sample,7:12);
rv_DRO2_LiAISON = DRO2_LiAISON(data_sample,7:12);
rv_LEO_real = LEO_real(data_sample,7:12);

% J2000转移至旋转系
[rv_DRO1_real_MCR,a_DRO1_real_MCR] = ECJ2k2Rot(t_sample,rv_DRO1_real,aux);
[rv_DRO1_LiAISON_MCR,a_DRO1_LiAISON_MCR] = ECJ2k2Rot(t_sample,rv_DRO1_LiAISON,aux);
[rv_DRO2_real_MCR,a_DRO2_real_MCR] = ECJ2k2Rot(t_sample,rv_DRO2_real,aux);
[rv_DRO2_LiAISON_MCR,a_DRO2_LiAISON_MCR] = ECJ2k2Rot(t_sample,rv_DRO2_LiAISON,aux);
[rv_LEO_real_MCR,a_LEO_real_MCR] = ECJ2k2Rot(t_sample,rv_LEO_real,aux);

% 定轨误差转移至轨道系
rvTC_DRO1_MCR = rv_DRO1_LiAISON_MCR-rv_DRO1_real_MCR;
rvTC_DRO1_TCR = T_TCR2TCO_eph(rvTC_DRO1_MCR,rv_DRO1_real_MCR,a_DRO1_real_MCR,FrameFlag);
rvTC_DRO2_MCR = rv_DRO2_LiAISON_MCR-rv_DRO2_real_MCR;
rvTC_DRO2_TCR = T_TCR2TCO_eph(rvTC_DRO2_MCR,rv_DRO2_real_MCR,a_DRO2_real_MCR,FrameFlag);

% % 星间相对误差
% rvTC_LiAISON_MCR = rv_DRO2_LiAISON_MCR-rv_DRO1_LiAISON_MCR;
% rvTC_LiAISON_TCR = T_TCR2TCO_eph(rvTC_LiAISON_MCR,rv_DRO1_LiAISON_MCR,a_DRO1_LiAISON_MCR,FrameFlag);
% rvTC_real_MCR = rv_DRO2_real_MCR-rv_DRO1_real_MCR;
% rvTC_real_TCR = T_TCR2TCO_eph(rvTC_real_MCR,rv_DRO1_real_MCR,a_DRO1_real_MCR,FrameFlag);
% error = rvTC_LiAISON_TCR-rvTC_real_TCR;

%% 画图
figure(1)
plot(t_sample/86400,rvTC_DRO1_TCR(:,1:3),'LineWidth',1.5)
title(['OD Pos Accuracy (',FrameFlag,', DRO1)'])
xlabel('t [day]'); ylabel('Position Accurcy [km]'); 
ylim([-0.03,0.03])
grid on; grid minor
legend('r_x','r_y','r_z')
set(gca,'FontSize',15,'fontname','times new roman');
figure(2)
plot(t_sample/86400,rvTC_DRO1_TCR(:,4:6),'LineWidth',1.5)
title(['OD Vel Accuracy (',FrameFlag,', DRO1)'])
xlabel('t [day]'); ylabel('Velocity Accurcy [km/s]'); 
ylim([-2,2]*1e-7)
grid on; grid minor
legend('v_x','v_y','v_z')
set(gca,'FontSize',15,'fontname','times new roman');

figure(3)
plot(t_sample/86400,rvTC_DRO2_TCR(:,1:3),'LineWidth',1.5)
title(['OD Pos Accuracy (',FrameFlag,', DRO2)'])
xlabel('t [day]'); ylabel('Position Accurcy [km]'); 
ylim([-0.03,0.03])
grid on; grid minor
legend('r_x','r_y','r_z')
set(gca,'FontSize',15,'fontname','times new roman');
figure(4)
plot(t_sample/86400,rvTC_DRO2_TCR(:,4:6),'LineWidth',1.5)
title(['OD Vel Accuracy (',FrameFlag,', DRO2)'])
xlabel('t [day]'); ylabel('Velocity Accurcy [km/s]'); 
ylim([-2,2]*1e-7)
grid on; grid minor
legend('v_x','v_y','v_z')
set(gca,'FontSize',15,'fontname','times new roman');

figure(5)
plot(t_sample/86400,rvTC_DRO1_TCR(:,1:3)-rvTC_DRO2_TCR(:,1:3),'LineWidth',1.5)
title(['OD Pos Accuracy (',FrameFlag,', Rel DRO)'])
xlabel('t [day]'); ylabel('Velocity Accurcy [km/s]'); 
ylim([-0.03,0.03])
grid on; grid minor
legend('r_x','r_y','r_z')
set(gca,'FontSize',15,'fontname','times new roman');

figure(6)
plot(t_sample/86400,rvTC_DRO1_TCR(:,4:6)-rvTC_DRO2_TCR(:,4:6),'LineWidth',1.5)
title(['OD Vel Accuracy (',FrameFlag,', Rel DRO)'])
xlabel('t [day]'); ylabel('Position Accurcy [km/s]'); 
ylim([-2,2]*1e-7)
grid on; grid minor
legend('v_x','v_y','v_z')
set(gca,'FontSize',15,'fontname','times new roman');
%% 判断theta与定轨精度的关系
% figure(100)
% r_LEODRO1_MCR = rv_DRO1_real_MCR(ceil(size_Sample/3):end,1:3)-rv_LEO_real_MCR(ceil(size_Sample/3):end,1:3);
% r_LEODRO1_TCR = T_TCR2TCO_eph(r_LEODRO1_MCR,rv_DRO1_real_MCR,a_DRO1_real_MCR,FrameFlag);
% theta_x = acos(r_LEODRO1_TCR(:,1)./sqrt(sum(r_LEODRO1_TCR.^2,2)));
% theta_y = acos(r_LEODRO1_TCR(:,2)./sqrt(sum(r_LEODRO1_TCR.^2,2)));
% % theta_z = acos(r_EarthDRO1_TCR(:,3)./sqrt(sum(r_EarthDRO1_TCR.^2,2)));
% % theta = atan2(rv_DRO1_real_MCR(ceil(size_Sample/3):end,2),rv_DRO1_real_MCR(ceil(size_Sample/3):end,1));
% % theta = acos(sum(r_DRO1.*r_EarthDRO1,2)./sqrt(sum(r_DRO1.^2,2))./sqrt(sum(r_EarthDRO1.^2,2)));
% 
% subplot(2,1,1)
% plot(theta_x, rvTC_DRO1_TCR(ceil(size_Sample/3):end,1),'.')
% title('OD Accuracy (TCR, DRO1)')
% xlabel('\theta_x [rad]'); ylabel('x [km]'); 
% xlim([0,pi])
% xticks([0, pi/4,pi/2,3*pi/4,pi]); xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'});
% grid on; grid minor
% set(gca,'FontSize',15,'fontname','times new roman');
% 
% subplot(2,1,2)
% plot(theta_y, rvTC_DRO1_TCR(ceil(size_Sample/3):end,2),'.')
% % title('Orbit Determination Accuracy (TCR)')
% xlabel('\theta_y [rad]'); ylabel('y [km]'); 
% xlim([0,pi]);
% xticks([0, pi/4,pi/2,3*pi/4,pi]); xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'});
% grid on; grid minor
% set(gca,'FontSize',15,'fontname','times new roman');
% 
% % subplot(2,1,2)
% % plot(theta_z, rvTC_DRO1_TCR(ceil(size_Sample/3):end,3),'.')
% % % title('Orbit Determination Accuracy (TCR)')
% % xlabel('\theta [rad]'); ylabel('z [km]'); 
% % xlim([0,pi]);
% % xticks([0, pi/4,pi/2,3*pi/4,pi]); xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'});
% % grid on; grid minor
% % set(gca,'FontSize',15,'fontname','times new roman');
