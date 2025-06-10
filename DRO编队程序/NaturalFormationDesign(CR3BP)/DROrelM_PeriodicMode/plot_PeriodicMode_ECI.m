% 画ECI中的分段编队与转移轨迹
clear
addpath('../../subF_eom(CR3BP)')
load('FloquetEig12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

%% 计算ECI下的相对运动
dt = 2*para.T0; % 积分时间                  
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;
y0_all = [100,150,200];
fig = figure(1);
hold off
for ii_index = 1:length(y0_all)

    x0_REL = 1/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';
    % 标称轨道及相对运动 的积分
    sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
    sol_sample = deval(sol,t_sample);
    xx_MCR_target = sol_sample(1:6,:)';
    xx_MCRLVLH_rel = sol_sample(7:12,:)';
    
    [xx_ECILVLH_rel,xx_ECI_target] = T_MCO2ECO_CRTBP(xx_MCRLVLH_rel,xx_MCR_target,t_sample,'LVLH',con.mu);
    p1 = plot(xx_ECILVLH_rel(:,1)*y0_all(ii_index)*con.r_norma,xx_ECILVLH_rel(:,2)*y0_all(ii_index)*con.r_norma,'LineWidth',1.5,'Color',[0, 114, 189]/255);
    hold on
    
    % 计算变轨轨迹
    ii_sample1 = 900;
    ii_sample2 = 1100;
    x0_MCR_target_trans = xx_MCR_target(ii_sample1,:);
    if ii_index == 1
        r_sep = y0_all(ii_index)/3*xx_MCRLVLH_rel(ii_sample1,1:3);
        r_sep_ECLLVLH_km = y0_all(ii_index)/3*xx_ECILVLH_rel(ii_sample1,1:3)*con.r_norma;
        r0f_MCRLVLH_rel_trans_km = [r_sep;y0_all(ii_index)*xx_MCRLVLH_rel(ii_sample2,1:3)]*con.r_norma;
        p01 = plot([0,r_sep_ECLLVLH_km(1)],[0,r_sep_ECLLVLH_km(2)],'LineWidth',1.5,'Color',[133, 133, 133]/255);
    else
        r0f_MCRLVLH_rel_trans_km = [y0_all(ii_index-1)*xx_MCRLVLH_rel(ii_sample1,1:3);y0_all(ii_index)*xx_MCRLVLH_rel(ii_sample2,1:3)]*con.r_norma;
    end
    t0 = t_sample(ii_sample1); tf = t_sample(ii_sample2); 
    dt_trans = (tf-t0);
    [~,~,~,xx_MCRLVLH_rel_trans_km,xx_MCR_target_trans] = forcedRelMotion(x0_MCR_target_trans,r0f_MCRLVLH_rel_trans_km,dt_trans,'LVLH',0,0);
    t_sample_trans = linspace(t0,tf,length(xx_MCRLVLH_rel_trans_km));
    [xx_ECILVLH_rel_trans_km,xx_ECI_target_trans] = T_MCO2ECO_CRTBP(xx_MCRLVLH_rel_trans_km,xx_MCR_target_trans,t_sample_trans,'LVLH',con.mu);

    p2 = plot(xx_ECILVLH_rel_trans_km(:,1),xx_ECILVLH_rel_trans_km(:,2),'LineWidth',1.5,'Color','r');
    p1f = plot(xx_ECILVLH_rel_trans_km(1,1),xx_ECILVLH_rel_trans_km(1,2),'rv');
    p10 = plot(xx_ECILVLH_rel_trans_km(end,1),xx_ECILVLH_rel_trans_km(end,2),'g^');
end


fig.Renderer = 'painters';
p0 = plot(0,0,'ks','MarkerSize',5);
bound = 3*max(y0_all);
% plot(bound*[-1,1],bound*[-cosd(60),cosd(60)],'Color',[1,1,1]*0.5);
% plot(bound*[1,-1],bound*[-cosd(60),cosd(60)],'Color',[1,1,1]*0.5)
poly1 = polyshape(bound*[1,-1,0],bound*[cosd(60),cosd(60),0]);
poly2 = polyshape(bound*[-1,1,0],bound*[-cosd(60),-cosd(60),0]);
py21 = plot(poly1,'LineWidth',1.5,'EdgeColor',[1,1,1]*0.5);
py21.FaceColor = [77, 190, 238]/255;
py21.FaceAlpha = 0.35;
py22 = plot(poly2,'LineWidth',1.5,'EdgeColor',[1,1,1]*0.5);
py22.FaceColor = [77, 190, 238]/255;
py22.FaceAlpha = 0.35;
xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); 
% title({'周期相对运动与主星可观测范围','(J2000 LVLH)'}); 

set(gca,'FontSize',15);

axis equal;
ylim(1.1*max(y0_all)*[-1,1])
xlim(1.5*max(y0_all)*[-1,1])
grid on; 
% grid minor
% legend([p0,p1,p21],{'A星','CRTBP参考轨道','A星可观测范围'},'Location','northeastoutside')
legend([p0,py21,p01,p2,p1,p10,p1f],{'A星','A星可观测范围','分离轨道','转移轨道','任务轨道','任务轨道起始点({\it\alpha}_{0})','任务轨道末端点({\it\alpha}_{\itf})'},'Location','northeastoutside')

set(gcf,'Color',[255,255,255]/255);
% export_fig PeModeEC.png -r600
