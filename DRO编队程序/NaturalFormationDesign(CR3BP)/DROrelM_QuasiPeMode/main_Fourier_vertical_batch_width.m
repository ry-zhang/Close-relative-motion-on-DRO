
%% CR3BP下DRO的相对运动法向拟周期解的不变环
% 验证法向拟周期轨迹的极值是否可准确落到theta0不变环上
% 结果表明，不能准确落到theta0不变环上，且DRO周期越大，误差越大。
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

flag = '4-2';
load(['generalSolFFT_all',flag,'.mat'])

opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
num = length(coe);
z_width_all = zeros(num,1);
diff_max_theta2_rel_all = zeros(num,1);
plot_flag = 0;
num_ii = 64;
theta2_0_all = linspace(0,2*pi,num_ii);
isplot = 1;

% parfor i_index = 1:num
for i_index = [1,100,200,300]
    disp(num2str(i_index))
    length_t = 2001;
    theta0_all_fit = linspace(0,2*pi,length_t);
    e5_hat_refit = real(iDFTmatrix_theta(coe(i_index).N,theta0_all_fit) * coe(i_index).c1_e5hat)';
    e6_hat_refit = real(iDFTmatrix_theta(coe(i_index).N,theta0_all_fit) * coe(i_index).c1_e6hat)';
    
    max_InvCircle_theta0 = zeros(num_ii,1);
    min_InvCircle_theta0 = zeros(num_ii,1);
    for ii_loop = 1:length(theta2_0_all)
        theta2_temp = theta2_0_all(ii_loop);
        e5_InvCircle = cos(theta2_temp)*e5_hat_refit + sin(theta2_temp)*e6_hat_refit;
        max_InvCircle_theta0(ii_loop) = max(e5_InvCircle(3,:));
        min_InvCircle_theta0(ii_loop) = min(e5_InvCircle(3,:));
    end
    maxmax_InvCir = max(max_InvCircle_theta0);
    minmax_InvCir = min(max_InvCircle_theta0);
    
    z_width_all(i_index) = (maxmax_InvCir-minmax_InvCir)/maxmax_InvCir;
    
    if isplot == 1
        length_t = 2001;
        theta0_all_fit = linspace(0,2*pi,length_t);

        % long-time evolution
        theta2_0 = pi/3;
        t_sample = linspace(0,150*para(i_index).T0,200000);
        e5_hat_refit_prop = real(iDFTmatrix_theta(coe(i_index).N,t_sample*2*pi/para(i_index).T0) * coe(i_index).c1_e5hat)';
        e6_hat_refit_prop = real(iDFTmatrix_theta(coe(i_index).N,t_sample*2*pi/para(i_index).T0) * coe(i_index).c1_e6hat)';
        theta2_temp = theta2_0+t_sample*2*pi/para(i_index).T2;
        rv_rel_normal = (cos(theta2_temp).*e5_hat_refit_prop + sin(theta2_temp).*e6_hat_refit_prop);
        ratio = 0.4/mean(abs(rv_rel_normal(3,:)));
        rv_rel_normal = ratio*rv_rel_normal;

        figure(3);
        p1 = plot(t_sample/para(i_index).T0,rv_rel_normal(3,:),'LineWidth',1.5); hold on

        upper_bound = ratio*interp1(theta2_0_all,max_InvCircle_theta0,mod(theta2_temp,2*pi));
        lower_bound = ratio*interp1(theta2_0_all,min_InvCircle_theta0,mod(theta2_temp,2*pi));
        pboundupper = plot(t_sample/para(i_index).T0,upper_bound,'LineWidth',1.5,'Color',[217, 83, 25]/255);
        plot(t_sample/para(i_index).T0,lower_bound,'LineWidth',1.5,'Color',[217, 83, 25]/255);
        
%         legend([p1,pboundupper],{'法向拟周期模态','上下界'},'Location','southoutside','NumColumns',2,'FontSize',17)
        grid on; box on; hold off
        ylim([-1,1])
        title(['\itT_{\rm0}\rm=',num2str(para(i_index).T0*con.T_norma_day,'%.2f'),'天'])
        set(gca,'FontSize',17)
        xlabel('\itt \rm[\itT_{\rm0}\rm]'); ylabel('\itz_L');
        exportgraphics(gcf,'normalMode.jpg','Resolution',600)
    end
end

%% 画图
figure(2)
plot([para.T0]*con.T_norma_day,z_width_all,'LineWidth',1.5);
xlabel('\itT_{\rm0} \rm[day]'); ylabel('\itz_{\rmwidth}');
set(gca,'FontSize',15); xlim([0,27])
grid on;

exportgraphics(gcf,'NMode_all.jpg','Resolution',600)