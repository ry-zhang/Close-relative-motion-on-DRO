% 计算不变环的极值波动

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

flag = '4-2';
load(['generalSolFFT_all',flag,'.mat'])

num = length(coe);
diff_max_theta2_abs_all = zeros(num,1);
diff_max_theta2_rel_all = zeros(num,1);
plot_flag = 1;
halfnum = coe(1).N/2;
bj = 1:halfnum;
aj = 0:halfnum;

colormap(zeros(num,3)); colormap parula;
color_all = colormap;
% 拟合颜色数据
T0_unif = linspace(para(num).T0,para(1).T0,num)';
fitType = fittype( 'cubicspline' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
fitobj_r = fit( T0_unif, color_all(:,1), fitType, opts );
fitobj_g = fit( T0_unif, color_all(:,2), fitType, opts );
fitobj_b = fit( T0_unif, color_all(:,3), fitType, opts );
ratio_all = zeros(100,1);
figure(100); clf;
for i_index = 1:num
    T0 = para(i_index).T0;
    color_T0 = [feval(fitobj_r,T0),feval(fitobj_g,T0),feval(fitobj_b,T0)];
    subplot(1,2,1)
    ratio = 1/imag(coe(i_index).c1_e5hat(halfnum,3));
    ratio_all(i_index) = ratio;
    semilogy(bj,ratio*imag(coe(i_index).c1_e5hat(halfnum:-1:1,3)),'.','MarkerSize',10,...
        'Color',color_T0);
    hold on
    grid on; grid minor; set(gca,'FontSize',15)
    xlim([0,25]); ylim([1e-15,1.5e0]);
    xlabel('\itj'); 
    ylabel('$b_j$ of $\hat{e}_5$','interpreter','latex');
    
    subplot(1,2,2)
    semilogy(aj,abs(ratio*real(coe(i_index).c1_e6hat(halfnum+1:-1:1,3))),'.','MarkerSize',10,...
        'Color',color_T0);
    hold on
    grid on; grid minor; set(gca,'FontSize',15)
    xlim([0,25]); ylim([1e-15,1.5e0]);
    xlabel('\itj'); 
    ylabel('$a_j$ of $\hat{e}_6$','interpreter','latex');
end
% b = colorbar;
% % b.YTick = [1,5,10,15,20];
% % b.YTickLabel = ;
% b.Limits = [0.1,27];
% b.Label.String = '\itT_{\rm0} \rm[day]';
% b.Label.FontSize = 15;
% caxis(b.Limits); % 更改颜色值的上下限

% set(gcf,'Color',[255,255,255]/255);
exportgraphics(gcf,'coeff.jpg','Resolution',600)

%% 画法向拟周期振动轨迹
for i_index = [1,100,200,300]
    disp(num2str(i_index))
    length_t = 2001;
    theta0_all_fit = linspace(0,2*pi,length_t);
    e5_hat_refit = real(iDFTmatrix_theta(coe(i_index).N,theta0_all_fit) * coe(i_index).c1_e5hat)';
    e6_hat_refit = real(iDFTmatrix_theta(coe(i_index).N,theta0_all_fit) * coe(i_index).c1_e6hat)';
    
    % long-time evolution
    theta2_0 = pi/3;
    t_sample = linspace(0,150*para(i_index).T0,200000);
    e5_hat_refit_prop = real(iDFTmatrix_theta(coe(i_index).N,t_sample*2*pi/para(i_index).T0) * coe(i_index).c1_e5hat)';
    e6_hat_refit_prop = real(iDFTmatrix_theta(coe(i_index).N,t_sample*2*pi/para(i_index).T0) * coe(i_index).c1_e6hat)';
    theta2_temp = theta2_0+t_sample*2*pi/para(i_index).T2;
    rv_rel_normal = (cos(theta2_temp).*e5_hat_refit_prop + sin(theta2_temp).*e6_hat_refit_prop);
    rv_rel_normal = 0.4/mean(abs(rv_rel_normal(3,:)))*rv_rel_normal;

    figure(3);
    p1 = plot(t_sample/para(i_index).T0,rv_rel_normal(3,:),'LineWidth',1.5);
    grid on; box on
    ylim([-1,1])
    title(['\itT_{\rm0}\rm=',num2str(para(i_index).T0*con.T_norma_day,'%.2f'),'天'])
    set(gca,'FontSize',17)
    xlabel('\itt \rm[\itT_{\rm0}\rm]'); ylabel('\itz_L\k_{\rm2}');
    
    exportgraphics(gcf,'normalMode.jpg','Resolution',600)
end
