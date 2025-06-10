clear
addpath('../../subF_eom(CR3BP)')
load('FloquetEig12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

%% y0
dt = 1*para.T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;
y0_all = [100,50,20,-20,-50,-100];
% y0_all = [200,90,20,-20,-90,-200];
% y0_all = 130;
% y0_all = [-5,-2.5,-1,1,2.5,5];
fig = figure(1);
hold off
for ii_index = 1:length(y0_all)
% for ii_index = 3
    x0_REL = y0_all(ii_index)/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';
    % 标称轨道及相对运动 的积分
    sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
    sol_sample = deval(sol,t_sample);
    abs_motion_M = [sol_sample(1:3,:);sol_sample(4:6,:)];
    rel_motion_L_linear = sol_sample(7:12,:);
    rel_motion_M_linear = T_TCO2TCR_CR3BP(rel_motion_L_linear',abs_motion_M','LVLH',con.mu)';
    
    plot(rel_motion_L_linear(1,:)*con.r_norma,rel_motion_L_linear(2,:)*con.r_norma,'LineWidth',1.5);
%     plot(rel_motion_M_linear(1,:)*con.r_norma,rel_motion_M_linear(2,:)*con.r_norma,'LineWidth',1.5);
%     rel_motion_VVLH_linear = T_TCR2TCO_CR3BP(rel_motion_M_linear',abs_motion_M','VVLH',con.mu)';
%     plot(rel_motion_VVLH_linear(1,:)*con.r_norma,rel_motion_VVLH_linear(3,:)*con.r_norma,'LineWidth',1.5);

    hold on
end
fig.Renderer = 'painters';
plot(0,0,'ks','MarkerSize',5)

plot(1.2*[min(y0_all),max(y0_all)]/2.775,1.2*[min(y0_all),max(y0_all)],'k--')
plot(1.2*[min(y0_all),max(y0_all)]/2.775,1.2*[max(y0_all),min(y0_all)],'k--')
xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); title('周期相对运动 (MCR LVLH)'); 
% text(0.4,0.1,'\leftarrow Leader','FontSize',15)
% plot(1.2*[min(y0_all),max(y0_all)],1.2*[min(y0_all),max(y0_all)]/2.775,'k--')
% plot(1.2*[max(y0_all),min(y0_all)],1.2*[min(y0_all),max(y0_all)]/2.775,'k--')
% xlabel('x [km]'); ylabel('y [km]'); title('Periodic Relative Orbit (VVLH frame)'); 
% text(0,-0.8,{'\uparrow','Leader'},'FontSize',15,'HorizontalAlignment','Center')
set(gca,'FontSize',15);

axis equal;
grid on; grid minor
str = strings(1,length(y0_all));
for ii = 1:length(y0_all)
    str(ii) = ['scale = ',num2str(y0_all(ii)),' km'];
%     str(ii) = ['y(0)=',num2str(y0_all(ii)),'km'];
end
% legend(str,'Location','north')
legend(str,'Location','northeastoutside')
% title('巡航卫星(L frame)')
set(gcf,'Color',[255,255,255]/255);
export_fig PeModeMul.png -r600

%% k0
dt = 1*para.T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;
% y0_all = [300,150,50,-50,-150,-300];
% y0_all = 130;
% y0_all = [-5,-2.5,-1,1,2.5,5];
k0_all = [1,0.5,0.2,-0.2,-0.5,-1]*1e-3;
% k0_all = -0.9*1e-5;
% y0_all = k0_all;
fig = figure(1);
hold off
for ii_index = 1:length(k0_all)
%     x0_REL = y0_all(ii_index)/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';
    x0_REL = k0_all(ii_index)*Sol_linear.vec3';
    % 标称轨道及相对运动 的积分
    sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
    sol_sample = deval(sol,t_sample);
    abs_motion_M = [sol_sample(1:3,:);sol_sample(4:6,:)];
    rel_motion_L_linear = sol_sample(7:12,:);
    rel_motion_M_linear = T_TCO2TCR_CR3BP(rel_motion_L_linear',abs_motion_M','LVLH',con.mu)';
    
    plot(rel_motion_L_linear(1,:)*con.r_norma,rel_motion_L_linear(2,:)*con.r_norma,'LineWidth',1.5);
%     rel_motion_VVLH_linear = T_TCR2TCO_CR3BP(rel_motion_M_linear',abs_motion_M','VVLH',con.mu)';
%     plot(rel_motion_VVLH_linear(1,:)*con.r_norma,rel_motion_VVLH_linear(3,:)*con.r_norma,'LineWidth',1.5);

    hold on
end
fig.Renderer = 'painters';
plot(0,0,'ks','MarkerSize',5)
ratio = con.r_norma*0.6;
% ratio = 1.1;
plot(ratio*[min(k0_all),max(k0_all)]/3.05,ratio*[min(k0_all),max(k0_all)],'k--')
plot(ratio*[min(k0_all),max(k0_all)]/3.05,ratio*[max(k0_all),min(k0_all)],'k--')
xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); 
title('周期相对运动 (坐标系\itL\rm)'); 
text(8,3,'\leftarrow 主星','FontSize',15)
% title('periodic mode (frame L)'); 
% text(8,3,'\leftarrow chief','FontSize',15)
set(gca,'FontSize',13);

axis equal;
ylim(ratio*[min(k0_all),max(k0_all)])
xlim([-200,200])
grid on; grid minor
str = strings(1,length(k0_all));
% for ii = 1:length(k0_all)
%     str(ii) = ['x(0)=',num2str(k0_all(ii)),'km'];
% %     str(ii) = ['y(0)=',num2str(y0_all(ii)),'km'];
% end
for ii = 1:length(k0_all)
%     str(ii) = ['x(0)=',num2str(y0_all(ii)),'km'];
    str(ii) = ['\itk_{\rm0}\rm=',num2str(k0_all(ii))];
%     str(ii) = ['y(0)=',num2str(y0_all(ii)),'km'];
end
% legend(str,'Location','north')
legend(str,'Location','northeastoutside')
% title('巡航卫星(L frame)')


set(gcf,'Color',[255,255,255]/255);
export_fig PeModeMul.png -r600