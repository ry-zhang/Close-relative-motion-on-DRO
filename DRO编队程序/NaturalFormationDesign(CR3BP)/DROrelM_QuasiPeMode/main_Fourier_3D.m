%% 三体轨道中的平面+法向拟周期解
% 2021-8-30
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
clear
addpath('../subF_eom(CR3BP)')

format longg
format compact

load('generalSolFFT_12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
length_t = 201;

%% 轨道积分分析
color1 = [224, 26, 255]/255;
color2 = [0, 255, 255]/255;
color_diff = (color1-color2)-sign(color1-color2)*1e-6;
num_ii = 17;
theta1_0_all = [0.625,0.625,0.625,0.625]*pi;
theta2_0_all = [0.25,0.875,1.25,1.75]*pi;
dt1 = 1*para.T0;
t_sample2 = linspace(0,dt1,length_t);

e1_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e1hat)';
e2_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e2hat)';
e5_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e5hat)';
e6_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample2*2*pi/para.T0) * coe.c1_e6hat)';
for ii_loop = 1:length(theta2_0_all)    
    theta1_0_temp = theta1_0_all(ii_loop);    
    theta2_0_temp = theta2_0_all(ii_loop);
    theta1_temp = mod(theta1_0_temp+t_sample2*2*pi/para.T1,2*pi);
    rv_rel_planar = cos(theta1_temp).*e1_hat_refit + sin(theta1_temp).*e2_hat_refit;
    theta2_temp = mod(theta2_0_temp+t_sample2*2*pi/para.T2,2*pi);
    rv_rel_vertical = cos(theta2_temp).*e5_hat_refit + sin(theta2_temp).*e6_hat_refit;
    rel_motion_temp = rv_rel_planar+rv_rel_vertical;
    
    
    figure(3);
    subplot(1,4,rem(ii_loop-1,4)+1)
    plot3(rel_motion_temp(1,:),rel_motion_temp(2,:),rel_motion_temp(3,:),'LineWidth',1.5); hold on
    plot3(0,0,0,'bo'); hold off
    grid on; axis equal; box on
    set(gcf,'Renderer','painters')
    xlim([-0.4,0.4]); ylim([-0.5,0.5]); zlim([-0.8,0.8]);
    set(gca,'FontSize',13)
    xlabel('t'); ylabel('z');
    title(['\theta_1(0)=',num2str(theta1_0_temp/pi),'\pi, ','\theta_2(0)=',num2str(theta2_0_temp/pi),'\pi'])
    pause(0.3)
end


