clear
addpath('../../subF_eom(CR3BP)')
load('FloquetEig12_sy.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

%% 计算ECI下的相对运动
dt = 2*para.T0; % 积分时间                  
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;
% y0_all = [100,50,20,-20,-50,-100];
% y0_all = [-5,-2.5,-1,1,2.5,5];
k0_all = [1,0.5,0.2,-0.2,-0.5,-1]*1e-3;
fig = figure(1);
hold off
for ii_index = 1:length(k0_all)
    x0_REL = k0_all(ii_index)*Sol_linear.vec3';
    % 标称轨道及相对运动 的积分
    sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
    sol_sample = deval(sol,t_sample);
    xx_MCR_target = sol_sample(1:6,:)';
    xx_MCRLVLH_rel = sol_sample(7:12,:)';
    
    [xx_ECILVLH_rel,xx_ECI_target] = T_MCO2ECO_CRTBP(xx_MCRLVLH_rel,xx_MCR_target,t_sample,'LVLH',con.mu);
    p1 = plot(xx_ECILVLH_rel(:,1)*con.r_norma,xx_ECILVLH_rel(:,2)*con.r_norma,'LineWidth',1.5);
    hold on
end
fig.Renderer = 'painters';
plot(0,0,'ks','MarkerSize',5)

% plot(1.2*[min(y0_all),max(y0_all)]/2.775,1.2*[min(y0_all),max(y0_all)],'k--')
% plot(1.2*[min(y0_all),max(y0_all)]/2.775,1.2*[max(y0_all),min(y0_all)],'k--')
xlabel('\itx_{LE} \rm[km]'); ylabel('\ity_{LE} \rm[km]'); title('平面周期模态'); 
% text(0.4,0.1,'\leftarrow Leader','FontSize',15)
% plot(1.2*[min(y0_all),max(y0_all)],1.2*[min(y0_all),max(y0_all)]/2.775,'k--')
% plot(1.2*[max(y0_all),min(y0_all)],1.2*[min(y0_all),max(y0_all)]/2.775,'k--')
% xlabel('x [km]'); ylabel('y [km]'); title('Periodic Relative Orbit (VVLH frame)'); 
% text(0,-0.8,{'\uparrow','Leader'},'FontSize',15,'HorizontalAlignment','Center')
set(gca,'FontSize',15);

axis equal;
grid on; grid minor
str = strings(1,length(k0_all));
for ii = 1:length(k0_all)
%     str(ii) = ['scale = ',num2str(y0_all(ii)),' km'];
    str(ii) = ['\itk_{\rm0}\rm=',num2str(k0_all(ii))];
end
% legend(str,'Location','north')
legend(str,'Location','northeastoutside')
% title('巡航卫星(L frame)')

set(gcf,'Color',[255,255,255]/255);
% export_fig PeModeMulEC.png -r600

%% 计算旋转坐标系DRO下的跟随编队在ECI下的不同,_f代表follower
dphase = para.T0/6;
t_sample_f = linspace(0+dphase,dt+dphase,length_t);
sol_f = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 dt+dphase], x0_DRO_M_3d, opts);
xx_MCR_chaser = deval(sol_f,t_sample_f);
xx_EMR_chaser = xx_MCR_chaser; xx_EMR_chaser(1:3,:) = xx_MCR_chaser(1:3,:)-[-1,0,0]';
xx_ECI_chaser = synodic2inertial(xx_EMR_chaser,t_sample)';

figure(2)
plot(xx_ECI_target(:,1),xx_ECI_target(:,2),'LineWidth',1.5); % 主星ECI下的轨道
hold on
plot(xx_ECI_chaser(:,1),xx_ECI_chaser(:,2),'LineWidth',1.5); % 副星ECI下的轨道
hold off
set(gca,'FontSize',15);
xlabel('\itx_{I} \rm[LU]'); ylabel('\ity_{I} \rm[LU]'); 
% title('DRO')
axis equal; grid on;
legend('主航天器','副航天器','Location','northeastoutside')
set(gcf,'Color',[255,255,255]/255);
export_fig ECI.png -r600
