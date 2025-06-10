%% CR3BP下DRO的相对运动法向拟周期解的不变环
% 2021-8-30
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../../subF_eom(CR3BP)')
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

format longg
format compact

% color1 = [224, 26, 255]/255;
% color2 = [0, 255, 255]/255;
% color1 = [0.2422, 0.1504, 0.6603];
% color2 = [0.9769, 0.9839, 0.0805];
% color_diff = (color1-color2)-sign(color1-color2)*1e-6;
num_ii = 64;
% map = color2+linspace(0,1,num_ii)'*color_diff;
colormap(zeros(num_ii,3));  colormap parula;
color_all = colormap;

load('generalSolFFT_12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

length_t = 2001;
theta0_all_fit = linspace(0,2*pi,length_t);
e5_hat_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e5hat)';
e6_hat_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e6hat)';

%% Invariant circles with constant \theta_2
% num_ii = 17;
figure(1); clf; hold on; 
theta2_0_all = linspace(0,2*pi,num_ii);
max_InvCircle_theta0 = zeros(num_ii,1);
min_InvCircle_theta0 = zeros(num_ii,1);
for ii_loop = 1:length(theta2_0_all)
    theta2_temp = theta2_0_all(ii_loop);
    e5_InvCircle = cos(theta2_temp)*e5_hat_refit + sin(theta2_temp)*e6_hat_refit;
    max_InvCircle_theta0(ii_loop) = max(e5_InvCircle(3,:));
    min_InvCircle_theta0(ii_loop) = min(e5_InvCircle(3,:));
%     plot(theta2_temp*ones(1,length_t),e5_InvCircle(3,:),'Color',color2+color_diff*ii_loop/num_ii,'LineWidth',1.5)
    plot(theta2_temp*ones(1,length_t),e5_InvCircle(3,:),'Color',color_all(ii_loop,:),'LineWidth',1.5)
    grid on; box on
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    set(gca,'FontSize',13)
    xlabel('\theta_2'); ylabel('\itz_L/k_{\rm2}');
    xlim([0,2*pi])
%     pause(0.3)
end
% title('Invariant curves of \theta_0')
title('不同的\theta_2对应的\theta_0不变环')
b = colorbar;
b.YTick = [0,pi/2,pi,3*pi/2,2*pi];
b.YTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
b.Limits = [0,2*pi];
b.Label.String = '\theta_2';
b.Label.FontSize = 15;
caxis([0 2*pi]); % 更改颜色值的上下限
% title('Invariant circles with constant \theta_1')
set(gca,'FontSize',15)
exportgraphics(gcf,'theta0.jpg','Resolution',600)

%% Invariant circles with constant \theta_0
figure(2); clf; hold on; 
theta0_all = linspace(0,2*pi,num_ii);
theta2_sam = linspace(0,2*pi,length_t);
% 记录不变环的极值
max_InvCircle_theta2 = zeros(num_ii,1);
min_InvCircle_theta2 = zeros(num_ii,1);
for ii_loop = 1:num_ii
    theta0_temp = theta0_all(ii_loop);
    
    e5_hat_temp = real(iDFTmatrix_theta(coe.N,theta0_temp) * coe.c1_e5hat)';
    e6_hat_temp = real(iDFTmatrix_theta(coe.N,theta0_temp) * coe.c1_e6hat)';
    e5_InvCircle = e5_hat_temp*cos(theta2_sam) + e6_hat_temp*sin(theta2_sam);
    
%     plot(theta0_temp*ones(1,length_t),e5_InvCircle(3,:),'Color',color2+color_diff*ii_loop/num_ii,'LineWidth',1.5)
    plot(theta0_temp*ones(1,length_t),e5_InvCircle(3,:),'Color',color_all(ii_loop,:),'LineWidth',1.5)
    max_InvCircle_theta2(ii_loop) = max(e5_InvCircle(3,:));
    min_InvCircle_theta2(ii_loop) = min(e5_InvCircle(3,:));
    grid on; box on
    xticks([0 pi/2 pi 3*pi/2 2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0,2*pi])
    set(gca,'FontSize',13)
    xlabel('\theta_0'); ylabel('\itz_L/k_{\rm2}');
    hold on; 
%     pause(0.3)
end
% title('Invariant curves of \theta_2')
title('不同的\theta_0对应的\theta_2不变环')
b = colorbar;

b.YTick = [0,pi/2,pi,3*pi/2,2*pi];
b.YTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
b.Limits = [0,2*pi];
b.Label.String = '\theta_0';
b.Label.FontSize = 15;
caxis([0 2*pi]); % 更改颜色值的上下限
% title('Invariant circles with constant \theta_1')
set(gca,'FontSize',15)
exportgraphics(gcf,'theta2.jpg','Resolution',600)

%% long-time evolution
bound_upper = max(max_InvCircle_theta0);
bound_lower = min(min_InvCircle_theta0);
bound_upper2 = min(max_InvCircle_theta0);
bound_lower2 = max(min_InvCircle_theta0);
figure(3);
theta2_0 = pi/3;
t_sample = linspace(0,20*para.T0,200000);
e5_hat_refit_prop = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e5hat)';
e6_hat_refit_prop = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e6hat)';
theta2_temp = theta2_0+t_sample*2*pi/para.T2;
rv_rel_normal = cos(theta2_temp).*e5_hat_refit_prop + sin(theta2_temp).*e6_hat_refit_prop;
p1 = plot(t_sample*con.T_norma_day,rv_rel_normal(3,:),'LineWidth',1.5); hold on;
% pboundupper = plot(t_sample*con.T_norma_day,bound_upper*ones(size(t_sample)),'--','LineWidth',1.5,'Color',[217, 83, 25]/255);
% plot(t_sample*con.T_norma_day,bound_lower*ones(size(t_sample)),'--','LineWidth',1.5,'Color',[217, 83, 25]/255);
% pboundlower = plot(t_sample*con.T_norma_day,bound_upper2*ones(size(t_sample)),'--','LineWidth',1.5,'Color',[252, 207, 48]/255);
% plot(t_sample*con.T_norma_day,bound_lower2*ones(size(t_sample)),'--','LineWidth',1.5,'Color',[252, 207, 48]/255);

upper_bound = interp1(theta2_0_all,max_InvCircle_theta0,mod(theta2_temp,2*pi));
lower_bound = interp1(theta2_0_all,min_InvCircle_theta0,mod(theta2_temp,2*pi));
pboundupper = plot(t_sample*con.T_norma_day,upper_bound,'--','LineWidth',1.5,'Color',[217, 83, 25]/255);
plot(t_sample*con.T_norma_day,lower_bound,'--','LineWidth',1.5,'Color',[217, 83, 25]/255);

hold off
legend([p1,pboundupper],{'法向拟周期模态','上下界'},...
    'Location','northeast')
grid on; box on
xlim([min(t_sample*con.T_norma_day),max(t_sample*con.T_norma_day)])
ylim([-1.2,1.2])
set(gca,'FontSize',15)
xlabel('\itt \rm[days]'); ylabel('\itz_L/k_{\rm2}');

% set(gcf,'Color',[255,255,255]/255);
% export_fig NBoundTrajL.png -r600
exportgraphics(gcf,'NBoundTrajL.jpg','Resolution',600)
%% 寻找极值点
% [peaks_max,peaks_index_max] = findpeaks(rv_rel_normal(3,:));
% [peaks_min,peaks_index_min] = findpeaks(-rv_rel_normal(3,:)); peaks_min = -peaks_min;
% figure(1)
% plot(mod(t_sample(peaks_index_max)*2*pi/para.T2+theta2_0,2*pi),peaks_max,'*');
% plot(mod(t_sample(peaks_index_min)*2*pi/para.T2+theta2_0,2*pi),peaks_min,'*')
% figure(2)
% plot(mod(t_sample(peaks_index_max)*2*pi/para.T0,2*pi),peaks_max,'*')
% plot(mod(t_sample(peaks_index_min)*2*pi/para.T0,2*pi),peaks_min,'*')
