%% CR3BP下DRO的相对运动平面拟周期解的不变环
% 2021-8-30
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
clear
addpath('../../subF_eom(CR3BP)')

format longg
format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 15);%坐标轴
set(0,'defaultTextFontSize', 15);%文字
set(0,'defaultLineLineWidth',1.5)

load('generalSolFFT_12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

length_t = 131;
theta0_all_fit = linspace(0,2*pi,length_t);
e1_hat_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e1hat)';
e2_hat_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e2hat)';

%% Invariant circles with constant \theta_1
num_ii = 9;

% color1 = [224, 26, 255]/255;
% color2 = [0, 255, 255]/255;
% color1 = [0.2422, 0.1504, 0.6603];
% color2 = [0.9769, 0.9839, 0.0805];
% color_diff = (color1-color2)-sign(color1-color2)*1e-6;
% map = color2+linspace(0,1,num_ii*4)'*color_diff;
colormap(zeros(num_ii*4,3)); colormap parula;
color_all = colormap;

figure(1); clf; hold on; 
% theta1_0_all = linspace(0,2*pi,num_ii); jj_num = 0;
theta1_0_all = linspace(0,pi/2,num_ii); title(['\theta_1=[0,\pi/2]']); jj_num = 0;
% theta1_0_all = linspace(pi/2,pi,num_ii); title(['\theta_1=[\pi/2,\pi]']); jj_num = 1;
% theta1_0_all = linspace(pi,3*pi/2,num_ii); title(['\theta_1=[\pi,3\pi/2]']); jj_num = 2;
% theta1_0_all = linspace(3*pi/2,2*pi,num_ii); title(['\theta_1=[3\pi/2,2\pi]']); jj_num = 3;
for ii_loop = 1:length(theta1_0_all)
    theta1_0_temp = theta1_0_all(ii_loop);
    e1_InvCircle = cos(theta1_0_temp)*e1_hat_refit + sin(theta1_0_temp)*e2_hat_refit;
%     plot(e1_InvCircle(1,:),e1_InvCircle(2,:),'Color',color2+color_diff*(ii_loop+jj_num*(num_ii-1)+1)/(4*(num_ii-1)+1),'LineWidth',1.5)
    plot(e1_InvCircle(1,:),e1_InvCircle(2,:),'Color',color_all(ii_loop+jj_num*(num_ii-1),:),'LineWidth',1.5)

%     p1 = plot(e1_InvCircle(1,1:3:end),e1_InvCircle(2,1:3:end),'Color',color2+color_diff*(ii_loop-1)/(num_ii-1));
%     p1.Color(4) = 0.05;
    axis equal; grid on; box on
    xlim([-0.4,0.4]); ylim([-0.5,0.5])
    set(gca,'FontSize',15)
    xlabel('\itx'); ylabel('\ity');
    
%     title(['\theta_1=',num2str(theta1_temp/pi),'\pi'])
%     pause(0.3)
end
% map = [color2; color2+color_diff/4; color2+color_diff/2; color2+3*color_diff/4; color1];

b = colorbar;
b.YTick = [0,pi/2,pi,3*pi/2,2*pi];
b.YTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
b.Limits = [0,2*pi];
b.Label.String = '\theta_1';
b.Label.FontSize = 15;
caxis([0 2*pi]); % 更改颜色值的上下限
% title('Invariant circles with constant \theta_1')
set(gca,'FontSize',15)
xlabel('\itx_L/k_{\rm1}'); ylabel('\ity_L/k_{\rm1}');
exportgraphics(gcf,'InvariantCurveTheta1.jpg','Resolution',600)

%% Invariant circles with constant \theta_0
num_ii = 9;

color1 = [0.2422, 0.1504, 0.6603];
color2 = [0.9769, 0.9839, 0.0805];
color_diff = (color1-color2)-sign(color1-color2)*1e-6;
map = color2+linspace(0,1,num_ii*4)'*color_diff;
colormap(map); colormap parula;
color_all = colormap;

figure(2); clf; hold on; 
% theta0_all = linspace(0,2*pi,num_ii);
theta0_all = linspace(0,pi,num_ii); title(['\theta_0=[0,\pi]']); jj_num = 0;
% theta0_all = linspace(pi,2*pi,num_ii); title(['\theta_0=[\pi,2\pi]']); jj_num = 1;
theta1_sam = linspace(0,2*pi,length_t);

for ii_loop = 1:length(theta0_all)
    theta0_temp = theta0_all(ii_loop);
%     t_temp = theta0_temp/2/pi*para.T0;
    e1_hat_temp = real(iDFTmatrix_theta(coe.N,theta0_temp) * coe.c1_e1hat)';
    e2_hat_temp = real(iDFTmatrix_theta(coe.N,theta0_temp) * coe.c1_e2hat)';
    e1_InvCircle = e1_hat_temp*cos(theta1_sam) + e2_hat_temp*sin(theta1_sam);
%     plot(e1_InvCircle(1,:),e1_InvCircle(2,:),'Color',color2+color_diff*(ii_loop+jj_num*(num_ii-1)+1)/(2*(num_ii-1)+1),'LineWidth',1.5)
    plot(e1_InvCircle(1,:),e1_InvCircle(2,:),'Color',color_all(ii_loop*2+jj_num*(num_ii*2-1),:),'LineWidth',1.5)

    axis equal; grid on; box on
    xlim([-0.4,0.4]); ylim([-0.5,0.5])
    set(gca,'FontSize',15)
    xlabel('\itx_L/k_{\rm1}'); ylabel('\ity_L/k_{\rm1}');
    hold on; 
%     title(['\theta_0=',num2str(theta1_temp/pi),'\pi'])
%     pause(0.3)
end
b = colorbar;
b.YTick = [0,pi/2,pi,3*pi/2,2*pi];
b.YTickLabel = {'0','\pi/2','\pi','3\pi/2','2\pi'};
b.Limits = [0,2*pi];
b.Label.String = '\theta_0';
b.Label.FontSize = 15;
caxis([0 2*pi]); % 更改颜色值的上下限
% title('Invariant circles with constant \theta_0')
exportgraphics(gcf,'InvariantCurveTheta0.jpg','Resolution',600)


%% 渲染边界图
num_bound = 2;
Bound_in = [];
Bound_out = [];
theta1_0_in = [linspace(pi/2-pi/20,pi/2+pi/20,num_bound),linspace(3*pi/2-pi/20,3*pi/2+pi/20,num_bound)];
theta1_0_out = [linspace(-pi/20,pi/20,num_bound),linspace(pi-pi/20,pi+pi/20,num_bound)];
for ii_loop = 1:length(theta1_0_in)
    e1_InvCircle_in = cos(theta1_0_in(ii_loop))*e1_hat_refit + sin(theta1_0_in(ii_loop))*e2_hat_refit;
    Bound_in = [Bound_in,polyshape(e1_InvCircle_in(1,:),e1_InvCircle_in(2,:))];
    e1_InvCircle_out = cos(theta1_0_out(ii_loop))*e1_hat_refit + sin(theta1_0_out(ii_loop))*e2_hat_refit;
    Bound_out = [Bound_out,polyshape(e1_InvCircle_out(1,:),e1_InvCircle_out(2,:))];
    
    figure(4)
    plot(e1_InvCircle_in(1,:),e1_InvCircle_in(2,:),'Color',[166, 166, 166]/255); hold on;
    plot(e1_InvCircle_in(1,:),e1_InvCircle_in(2,:),'.','Color',0.3*[166, 166, 166]/255);
    pinner = plot(polyshape(e1_InvCircle_in(1,:),e1_InvCircle_in(2,:))); hold off;
    pinner.FaceColor = [0.5,0.5,0.5]; pinner.EdgeColor = 'none';
    axis equal; xlim([-0.401,0.401]); ylim([-0.5,0.5]);
    set(gca,'XTick',[-0.4,0,0.4])
    axis off
    export_fig innerBound.png -r600
    figure(5);
    plot(e1_InvCircle_out(1,:),e1_InvCircle_out(2,:),'Color',[166, 166, 166]/255);hold on;
    plot(e1_InvCircle_out(1,:),e1_InvCircle_out(2,:),'.','Color',0.3*[166, 166, 166]/255);
    pouter = plot(polyshape(e1_InvCircle_out(1,:),e1_InvCircle_out(2,:))); hold off;
    pouter.FaceColor = [0.5,0.5,0.5];pouter.EdgeColor = 'none';
    axis equal; xlim([-0.401,0.401]); ylim([-0.5,0.5]);
    set(gca,'XTick',[-0.4,0,0.4])
    axis off
    export_fig outerrBound.png -r600
end
polyin = intersect(Bound_in); % 求交集
polyout = union(Bound_out); % 并集
polyall = subtract(polyout,polyin); % 差集

figure(4);
% plot(polyin.Vertices(:,1),polyin.Vertices(:,2),'Color',); hold on;
pinner = plot(polyin); hold on;
pinner.FaceAlpha = 0.35; pinner.FaceColor = [0.5,0.5,0.5]; 
pinner.EdgeColor = [166, 166, 166]/255; pinner.LineWidth = 1.5;
plot(polyin.Vertices(:,1),polyin.Vertices(:,2),'.','Color',0.3*[166, 166, 166]/255); hold off
axis equal; xlim([-0.401,0.401]); ylim([-0.5,0.5]); axis off
% title('inner boundary'); xlabel('\itx_L'); ylabel('\ity_L');
set(gcf,'Color',[255,255,255]/255); grid on
export_fig innerBound.png -r600

figure(5);
pouter = plot(polyout); hold on;
pouter.FaceAlpha = 0.35; pouter.FaceColor = [0.5,0.5,0.5]; 
pouter.EdgeColor = [166, 166, 166]/255; pouter.LineWidth = 1.5;
plot(polyout.Vertices(:,1),polyout.Vertices(:,2),'.','Color',0.3*[166, 166, 166]/255); hold off
axis equal; xlim([-0.401,0.401]); ylim([-0.5,0.5]); axis off
% title('inner boundary'); xlabel('\itx_L'); ylabel('\ity_L');
set(gcf,'Color',[255,255,255]/255); grid on
export_fig outerrBound.png -r600

figure(5);
pouter = plot(polyall); hold on;
pouter.FaceAlpha = 0.35; pouter.FaceColor = [0.5,0.5,0.5]; 
pouter.EdgeColor = [166, 166, 166]/255; pouter.LineWidth = 1.5;
plot(polyall.Vertices(:,1),polyall.Vertices(:,2),'.','Color',0.3*[166, 166, 166]/255); hold off
axis equal; xlim([-0.401,0.401]); ylim([-0.5,0.5]); axis off
% title('inner boundary'); xlabel('\itx_L'); ylabel('\ity_L');
set(gcf,'Color',[255,255,255]/255); grid on
export_fig allBound.png -r600


figure(3);
pall = plot(polyall); 
pall.FaceAlpha = 0.7; pall.FaceColor = [0.5,0.5,0.5];
pall.EdgeColor = 'none';xlabel('\itx_L'); ylabel('\ity_L');
title('duration = 80\itT')
set(gcf,'Color',[255,255,255]/255);

theta1_0 = pi/3;
t_sample = linspace(0,7*para.T0,20000);
e1_hat_refit_prop = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e1hat)';
e2_hat_refit_prop = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e2hat)';
theta1_temp = theta1_0+t_sample*2*pi/para.T1;
rv_rel_planar = cos(theta1_temp).*e1_hat_refit_prop + sin(theta1_temp).*e2_hat_refit_prop;
hold on;
p1 = plot(rv_rel_planar(1,:),rv_rel_planar(2,:),'LineWidth',1.5,'Color',[1, 114, 189]/255);
% p1 = plot(rv_rel_planar(1,:),rv_rel_planar(2,:),'LineWidth',1.5);
p2 = plot(rv_rel_planar(1,1),rv_rel_planar(2,1),'g^');
p3 = plot(rv_rel_planar(1,end),rv_rel_planar(2,end),'rv');
hold off
% legend([pall,p1,p2,p3],{'Bounded area','Trajectory','Initial Position','Final Position'},...
%     'Location','northeastoutside')
legend([pall,p1,p2,p3],{'有界区域','轨迹','初始位置','末端位置'},...
    'Location','northeastoutside')
axis equal; grid on; box on
xlim([-0.4,0.4]); ylim([-0.5,0.5])
set(gca,'FontSize',15)
xlabel('\itx_L/k_{\rm1}'); ylabel('\ity_L/k_{\rm1}');

% title('duration = 80\itT')
set(gcf,'Color',[255,255,255]/255);
exportgraphics(gcf,'PBoundTrajL.jpg','Resolution',600)
