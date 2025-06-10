
%% CR3BP下DRO的相对运动平面拟周期解的不变环
% 2021-8-30
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
clear
warning off
addpath('../../subF_eom(CR3BP)')

format longg
format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

flag = '4-2';
load(['generalSolFFT_all',flag,'.mat'])

opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
num = length(coe);
xp_all = zeros(num,2);
yp_all = zeros(num,2);
plot_flag = 1;

% for i_index = 1:num
for i_index = 120
% for i_index = [1,100,200,300]
    disp(num2str(i_index))
    length_t = 201;
    theta0_all_fit = linspace(0,2*pi,length_t);
    e1_hat_refit = real(iDFTmatrix_theta(coe(i_index).N,theta0_all_fit) * coe(i_index).c1_e1hat)';
    e2_hat_refit = real(iDFTmatrix_theta(coe(i_index).N,theta0_all_fit) * coe(i_index).c1_e2hat)';

    %% 渲染边界图
    num_bound = 9;
    Bound_in = [];
    Bound_out = [];
    theta1_0_in = [linspace(pi/2-pi/8,pi/2+pi/8,num_bound),linspace(3*pi/2-pi/8,3*pi/2+pi/8,num_bound)];
    theta1_0_out = [linspace(-pi/8,pi/8,num_bound),linspace(pi-pi/8,pi+pi/8,num_bound)];
    for ii_loop = 1:length(theta1_0_in)
        e1_InvCircle_in = cos(theta1_0_in(ii_loop))*e1_hat_refit + sin(theta1_0_in(ii_loop))*e2_hat_refit;
        Bound_in = [Bound_in,polyshape(e1_InvCircle_in(1,:),e1_InvCircle_in(2,:))];
        e1_InvCircle_out = cos(theta1_0_out(ii_loop))*e1_hat_refit + sin(theta1_0_out(ii_loop))*e2_hat_refit;
        Bound_out = [Bound_out,polyshape(e1_InvCircle_out(1,:),e1_InvCircle_out(2,:))];
    end
    polyin = intersect(Bound_in); % 求交集
    polyout = union(Bound_out); % 并集
    polyall = subtract(polyout,polyin); % 差集
    % 计算边界的边界
    [~,yp1_index] = min(abs(polyin.Vertices(:,1)));
    [~,yp2_index] = min(abs(polyout.Vertices(:,1)));
    [~,xp1_index] = min(abs(polyin.Vertices(:,2)));
    [~,xp2_index] = min(abs(polyout.Vertices(:,2)));
    xp_all_i = abs([polyin.Vertices(xp1_index,1), polyout.Vertices(xp2_index,1)]);
    yp_all_i = abs([polyin.Vertices(yp1_index,2), polyout.Vertices(yp2_index,2)]);
    xp_all(i_index,:) = xp_all_i;
    yp_all(i_index,:) = yp_all_i;
    if plot_flag == 1
        figure(3);
        polyall_resize = polyall;
        factor_resize = 0.4/yp_all_i(2);
        polyall_resize.Vertices = factor_resize*polyall.Vertices;
        pall = plot(polyall_resize); 
        pall.FaceAlpha = 0.7; pall.FaceColor = [0.5,0.5,0.5];
        pall.EdgeColor = 'none';

        theta1_0 = pi/3;
        t_sample = linspace(0,6*para(i_index).T0,2000);
        e1_hat_refit_prop = real(iDFTmatrix_theta(coe(i_index).N,t_sample*2*pi/para(i_index).T0) * coe(i_index).c1_e1hat)';
        e2_hat_refit_prop = real(iDFTmatrix_theta(coe(i_index).N,t_sample*2*pi/para(i_index).T0) * coe(i_index).c1_e2hat)';
        theta1_temp = theta1_0+t_sample*2*pi/para(i_index).T1;
        rv_rel_planar = cos(theta1_temp).*e1_hat_refit_prop + sin(theta1_temp).*e2_hat_refit_prop;
        rv_rel_planar = factor_resize*rv_rel_planar;
        hold on;
        p1 = plot(rv_rel_planar(1,:),rv_rel_planar(2,:),'LineWidth',1.5,'Color',[1, 114, 189]/255);
%         p2 = plot(rv_rel_planar(1,1),rv_rel_planar(2,1),'g^');
%         p3 = plot(rv_rel_planar(1,end),rv_rel_planar(2,end),'rv');
        hold off
%         legend([pall,p1,p2,p3],{'Bounded area','Trajectory','Initial Position','Final Position'},...
%             'Location','northeastoutside')
        axis equal; grid on; box on
        xlim([-0.4001,0.4001]); ylim([-0.5001,0.5001])
        set(gca,'FontSize',17)
        xlabel('\itx_L/k_{\rm1}'); ylabel('\ity_L/k_{\rm1}');
        title(['\itT_{\rm0}\rm=',num2str(para(i_index).T0*con.T_norma_day,'%.2f'),'天'])

        pause(0.1)
%         set(gcf,'Color',[255,255,255]/255);
        exportgraphics(gcf,'fig.jpg','Resolution',600)
    end
end

%% 画有界区域的覆盖范围变化
x_width = (xp_all(:,2)-xp_all(:,1))./xp_all(:,2);
y_width = (yp_all(:,2)-yp_all(:,1))./yp_all(:,2);
figure(4)
plot([para.T0]*con.T_norma_day,x_width,[para.T0]*con.T_norma_day,y_width,'LineWidth',1.5)
legend('\itx_{\rmwidth}','\ity_{\rmwidth}','Location','north')
grid on; box on
set(gca,'FontSize',15)
xlabel('\itT_{\rm0}\rm [day]'); ylabel('拟周期模态几何参数');
xlim([0,27])
% set(gcf,'Color',[255,255,255]/255);
exportgraphics(gcf,'NoMode_all.jpg','Resolution',600)