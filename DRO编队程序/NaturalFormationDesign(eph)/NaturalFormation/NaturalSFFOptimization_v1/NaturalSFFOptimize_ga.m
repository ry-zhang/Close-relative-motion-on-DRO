clear
% close all
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg; format compact
warning off

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

load('FloquetEig12')
opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-12);

%% 常数与变量
% AnalyzeFlag = 1; % 1-分析,0-优化
% load('DRO2NaturForm90day_DE.mat');

AnalyzeFlag = 0; % 1-分析,0-优化

% 优化/计算区间
% 为设计转移轨道，可以把编队优化时间调长，这样编队时间段可以来回滑动
numDay = 100;
tspan_sec = [0,86400*numDay];
t_sample = linspace(tspan_sec(1),tspan_sec(2),1500);
t0_forma_all = 86400*100*(0:11);
% 参考轨道尺寸，1-对应于参考轨道星间距离[0.67,1] km
% scale_ref = 20; % 5-50km   近距离编队
scale_ref = 90; % 50-100km 中距离编队
% scale_ref = 1000; % >500km 远距离编队
iDRO = 3; 
FileName = ['NaturalSSF_DRO',num2str(iDRO),'_',num2str(numDay),'day_scale',num2str(scale_ref),'_DE'];

%% CR3BP下的DRO周期相对运动
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = -1/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';

% 周期轨道
tf = para.T0;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
length_t = 2000;
t_sample_CRTBP = linspace(0,tf,length_t);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
rvDRO_MCR_ref = sol_deval_CRTBP(:,1:6);

rRelMot_ref_L = sol_deval_CRTBP(:,7:9);
rRelMot_ref = rRelMot_ref_L(:,[2,3,1]);
rRelMot_ref(:,[2,3]) = -rRelMot_ref(:,[2,3]);

vRelMot_ref_L = sol_deval_CRTBP(:,10:12);
vRelMot_ref = vRelMot_ref_L(:,[2,3,1]);
vRelMot_ref(:,[2,3]) = -vRelMot_ref(:,[2,3]);

rRelMot_km_L = rRelMot_ref_L.*con.r_norma;

%% 拟合参考周期轨道
rRelMotion_km = rRelMot_ref.*con.r_norma;
vRelMotion_km = vRelMot_ref.*con.v_norma;

theta_all = mod(atan2(-rvDRO_MCR_ref(:,2),rvDRO_MCR_ref(:,1)),2*pi);
% theta_all(end) = theta_all(end);
% Set up fittype and options.
ft = fittype( 'fourier8' );
opts_fit = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_fit.Display = 'Off';
% opts_fit.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.00074919498648];

% Fit model to data.
[fitresult_x, gof_x] = fit( theta_all, rRelMotion_km(:,1), ft, opts_fit );
[fitresult_z, gof_z] = fit( theta_all, rRelMotion_km(:,3), ft, opts_fit );

[fitresult_vx, gof_vx] = fit( theta_all, vRelMotion_km(:,1), ft, opts_fit );
[fitresult_vz, gof_vz] = fit( theta_all, vRelMotion_km(:,3), ft, opts_fit );

%% ga优化星历下的自然编队构型
aux = []; % 加载星历、设置初始历元
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;

if AnalyzeFlag ~= 1
    % 星历积分初值
    
    % 太偏了，不用
    % aux.t0UTC  = datetime(2023, 1, 1, 0, 0, 0); % 初始历元
    % x0_MCR_target = [72687.2175909459 0 0 0 -0.540957166578942 0];
    switch iDRO
    case 1
        % DRO1
        load('DROephMultiRev2024.mat')
        aux.jd0 = jd0;
        x0_j2k_target = x0j2k';    
    case 2
        % DRO2
        aux.t0UTC  = datetime(2023, 1, 1, 0, 0, 0); % 初始历元
        x0_j2k_target = [380224.0673756608739495, 140818.4956395560002420, 42079.2399940182804130, ...
            -0.5874909506324951, 0.6787788399867067, 0.3426576490741002];% 星历积分初值（j2000）
    case 3
        % DRO3 2023.6-2027.1 <<1 hr
        aux.t0UTC = datetime('2023-06-14 07:36:12');
        x0_j2k_target = [260861.357166097,167404.871364748,78987.9106263715,...
            -0.797706314601184,0.975263733684134,0.530154969412835];
    case 4
        % DRO4 2023.6-2027.6  ~2.5 hr
        aux.t0UTC = datetime('2023-06-13 23:00:16');
        x0_j2k_target = [276577.520222802,144514.098050027,67426.0628386977,...
            -0.708023138517267,1.02870635679717,0.558158622414164];
    otherwise
        error('no exist DRO')
    end
    
else
    aux.t0UTC = t0UTC;
end

aux = initialize(aux); % 初始化
t0UTC = aux.t0UTC;
if exist('x0_MCR_target','var') && ~exist('x0_j2k_target','var')
    x0_j2k_target = T_Rot2ECJ2k(aux.jd0, x0_MCR_target,aux.C_Mat, 'MCEMR');
elseif ~exist('x0_MCR_target','var') && exist('x0_j2k_target','var')
    x0_MCR_target = T_ECJ2k2Rot(aux.jd0, x0_j2k_target, [0,0,0],aux.C_Mat, 'MCEMR');
end


% load FHL
% aux.jd0  = data.jd0; % 初始历元
% aux.t0UTC = [];
% x0_MCR_target = data.xMT_MCR; % 星历DRO初值

% 优化
if AnalyzeFlag == 1
    [xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);
    theta_DRO_eph = mod(atan2(-xx_MCR_target(:,2),xx_MCR_target(:,1)),2*pi);
else
    % 星历积分区间 
%     tspan_sec = [0,para.T0*con.T_ref_day*86400*1];% 
    % 主星星历积分
    [xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);
    
    % 计算参考轨道
    theta_DRO_eph = mod(atan2(-xx_MCR_target(:,2),xx_MCR_target(:,1)),2*pi);
    fit_x = scale_ref*feval(fitresult_x,theta_DRO_eph);
    fit_z = scale_ref*feval(fitresult_z,theta_DRO_eph);
    rRel_ref = [fit_x,fit_z];
    
    r_min = scale_ref*min(rRelMotion_km(:,[1,3])); r_max = scale_ref*max(rRelMotion_km(:,[1,3]));
    v_min = scale_ref*min(vRelMotion_km(:,[1,3])); v_max = scale_ref*max(vRelMotion_km(:,[1,3]));
    scale_factor = [1,1,1e6,1e6]; % 优化归一化参数
    r0 = [fit_x(1),fit_z(1)];
    v0 = scale_ref*[feval(fitresult_vx,theta_DRO_eph(1)),...
        feval(fitresult_vz,theta_DRO_eph(1))];
    lb = scale_factor.*[r0-(r_max-r_min)/4, v0-(v_max-v_min)/4];
    ub = scale_factor.*[r0+(r_max-r_min)/4, v0+(v_max-v_min)/4];
    
    fx = @(x)dist_DROformation(x,xx_MCR_target,a_MCR_target,...
            tspan_sec,t_sample,rRel_ref,scale_factor,aux);
    opts_DE = struct('Max_nfes',4*1000,'Display',1,...
            'Resize',0,'UseParallel',1);
%     options = optimoptions('ga','Display','iter','UseParallel',true,'MaxGenerations',100,'PopulationSize',63);

    num_opt = 1;
    x0_TC_optimize_all = zeros(num_opt,4);
    MaxDist_all = zeros(num_opt,1);
    disp(['   优化次数' , '    最大误差', '     所耗时间']);
    for ii_loop = 1:num_opt
        tic
        
        [fval,x0_TC_optim_scale] = PfhDE(fx,lb,ub,opts_DE);
        
%         [x0_TC_optimize,fval] = ga(@(x)dist_DROformation(x,...
%             xx_MCR_target,a_MCR_target,tspan_sec,t_sample,rRel_ref,scale_factor,aux),...
%             4,[],[],[],[],lb,ub,[],options);
        
        x0_TC_optimize_all(ii_loop,:) = x0_TC_optim_scale(end,:)./scale_factor;
        MaxDist_all(ii_loop) = fval(end);
        t_end = toc;
        disp(['      ',num2str(ii_loop), '      ',num2str(fval(end),'%07.4f'),...
            'km     ', num2str(t_end,'%04.2f'),'s']);
    end
    save(FileName,'MaxDist_all', 'x0_TC_optimize_all', 'tspan_sec',...
        't_sample', 't0UTC', 'x0_MCR_target', 'scale_ref')
end

%% 计算优化结果
[~,num] = min(MaxDist_all);
% num = 1;
scale_0 = 1; % 初值放缩
x0_TCO_rel_optimal = scale_0*x0_TC_optimize_all(num,:);
x0_TCR_chaser = [x0_TCO_rel_optimal(1),0,x0_TCO_rel_optimal(2),x0_TCO_rel_optimal(3),0,x0_TCO_rel_optimal(4)];
x0_MCR_chaser = T_TCO2TCR_eph(x0_TCR_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;
% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);
% 计算相对运动
rv_MCR_rel = xx_MCR_chaser-xx_MCR_target;
rv_TCO_rel = T_TCR2TCO_eph(rv_MCR_rel,xx_MCR_target,a_MCR_target,'LVLH');

%% LVLH坐标系下的相对运动
%---------月心旋转系下的LVLH-----------
figure(1)
plot(0,0,'ks','MarkerSize',5); hold on;
plot(scale_0*scale_ref*rRelMot_km_L(:,1), scale_0*scale_ref*rRelMot_km_L(:,2),'k','LineWidth',1.5); 
plot(rv_TCO_rel(:,1), rv_TCO_rel(:,2),'Color',[237, 177, 32]/255,'LineWidth',1.5);
plot(rv_TCO_rel(1,1), rv_TCO_rel(1,2),'g^');
plot(rv_TCO_rel(end,1), rv_TCO_rel(end,2),'rv');
legend('主星','CRTBP','星历','初始时刻','终端时刻','Location','northeastoutside');
box on; grid on; grid minor; hold off; 
axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]')
ylim(scale_ref*[-1.3,0.2]); xlim(scale_ref*[-1,1]);
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title('月心旋转系下的LVLH坐标系')

%---------地心惯性系下的LVLH-----------
figure(2)
% CRTBP的地心惯性系下的LVLH
rvDRO_ECR_CRTBP = [sol_deval_CRTBP(:,1:3)'-[-1,0,0]'; sol_deval_CRTBP(:,4:6)'];
rvDRO_ECI_CRTBP = synodic2inertial(rvDRO_ECR_CRTBP,t_sample_CRTBP)';
rRelMot_ECR_linear = T_TCO2TCR_CR3BP(scale_0*scale_ref*rRelMot_km_L,rvDRO_MCR_ref,'LVLH',con.mu)'; 
rRelMot_ECI = synodic2inertial(rRelMot_ECR_linear,t_sample_CRTBP)';
rv_TCO_ECI_f = T_TCI2TCO_E_CR3BP(rRelMot_ECI,rvDRO_ECI_CRTBP,t_sample_CRTBP,'LVLH',con.mu)';
% 星历的地心惯性系下的LVLH
sol_j2k_target = ode113(@(t, xx)eqom_geoMEMEJ2k(t, xx, aux), tspan_sec, x0_j2k_target , opts_ode);
[xx_j2k_target,a_j2k_target] = deval(sol_j2k_target,t_sample);
xx_j2k_target = xx_j2k_target'; a_j2k_target = a_j2k_target(1:3,:)';
% xx_j2k_target = T_Rot2ECJ2k(aux.jd0+t_sample/86400,xx_MCR_target,aux.C_Mat, 'MCEMR',0);
xx_j2k_chaser = T_Rot2ECJ2k(aux.jd0+t_sample/86400,xx_MCR_chaser,aux.C_Mat, 'MCEMR',0);
rv_RelMot_j2k = xx_j2k_chaser-xx_j2k_target;
rv_RelMot_j2k_LVLH = T_TCR2TCO_eph(rv_RelMot_j2k,xx_j2k_target,a_j2k_target,'LVLH');
rv_RelMot_j2k_LVLH = rv_RelMot_j2k_LVLH';

% 画图
plot(0,0,'ks','MarkerSize',5); hold on
plot(rv_TCO_ECI_f(1,:),rv_TCO_ECI_f(2,:),'k','LineWidth',1.5);
plot(rv_RelMot_j2k_LVLH(1,:),rv_RelMot_j2k_LVLH(2,:),'LineWidth',1.5,'Color',[237, 177, 32]/255);
plot(rv_RelMot_j2k_LVLH(1,1),rv_RelMot_j2k_LVLH(2,1),'g^');
plot(rv_RelMot_j2k_LVLH(1,end),rv_RelMot_j2k_LVLH(2,end),'rv');
hold off; 
box on; grid on; grid minor; 
axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]')
legend('主星','CRTBP','星历','初始时刻','终端时刻','Location','northeastoutside');
title('地心J2000下的LVLH坐标系')
set(gca,'FontSize',15);

%% DRO 绝对运动
%---------月心旋转系下的DRO-----------
f3 = figure(3);
p2 = plot3(xx_MCR_target(:,1),xx_MCR_target(:,2),xx_MCR_target(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5); hold on;
p3 = plot3(xx_MCR_target(1,1),xx_MCR_target(1,2),xx_MCR_target(1,3),'g^');
p4 = plot3(xx_MCR_target(end,1),xx_MCR_target(end,2),xx_MCR_target(end,3),'rv');
p1 = plot3(rvDRO_MCR_ref(:,1)*aux.LU,rvDRO_MCR_ref(:,2)*aux.LU,rvDRO_MCR_ref(:,3)*aux.LU,'k','LineWidth',1.5);
box on; grid on; grid minor; hold off;
axis equal;xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]'); zlabel('\itz \rm[km]')
% L = legend([p1 p2,p3,p4],'CRTBP','eph','初始时刻','终端时刻','Location','northeast');
L = legend([p1 p2,p3,p4],'CRTBP','eph','初始时刻','终端时刻','Location','northeastoutside');
% set(L,'box','off')
xlim([-1.2,1.2]*1e5)
set(gca,'FontSize',15);
title('DRO (MCR)')
% title('DRO (M Frame)')
set(gcf,'Renderer','painters')
view(0,90)

%---------月心旋转系下的DRO-----------
xx_j2k_target = T_Rot2ECJ2k(aux.jd0+t_sample/86400,xx_MCR_target,aux.C_Mat, 'MCEMR',0);
t_sample_temp = linspace(0,2*para.T0,2000);
dphase = -para.T0/4.8;
[~,rDRO_MCR_ref] = ode113(@(t,x)eom_abs3b(t,x,con.mu),t_sample_temp, x0_DRO_M_3d, opts_ode); % 注意这个旋转系与MCR不一样

% 星历下DRO与CRTBP下DRO在MCR坐标系下是相近的，但是在转化到ECI的时候，
% 星历下考虑的地月距离波动较大，而CRTBP下考虑的地月距离是恒定的，
% 因此二者在ECI下会差别过大，除非将CRTBP转化的时候采用实时的地月距离。
% 由于LVLH坐标系下是已经优化过的，因此在LVLH坐标系下的相对轨迹反而差别较小

rDRO_ECI_ref = synodic2inertial((rDRO_MCR_ref(:,1:3)'-[-1,0,0]'),dphase+t_sample_temp)'*con.r_norma;
f4 = figure(4);
set(f4,'name','星历DRO ECI')
% figure('color',[1 1 1],'name','星历DRO MCR')
p2 = plot3(xx_j2k_target(:,1),xx_j2k_target(:,2),xx_j2k_target(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5); hold on;
p3 = plot3(xx_j2k_target(1,1),xx_j2k_target(1,2),xx_j2k_target(1,3),'g^');
p4 = plot3(xx_j2k_target(end,1),xx_j2k_target(end,2),xx_j2k_target(end,3),'rv');
% p1 = plot3(rDRO_ECI_ref(:,1),rDRO_ECI_ref(:,2),rDRO_ECI_ref(:,3),'k','LineWidth',1.5);
box on; grid on; grid minor; hold off;
axis equal;xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]'); zlabel('\itz \rm[km]')
L = legend([p2,p3,p4],'CRTBP','eph','初始时刻','终端时刻','Location','northeastoutside');
% set(L,'box','off')
% xlim([-1.2,1.2]*1e5)
set(gca,'FontSize',15);
title('DRO (J2000)')
set(gcf,'Renderer','painters')
view(0,90)

%% 画拟合图
% fitresult_x_L = fit( theta_all, rRelMotion_km_L(:,1), ft, opts_fit );
% fitresult_y_L = fit( theta_all, rRelMotion_km_L(:,2), ft, opts_fit );
% figure(5)
% subplot(2,1,1)
% hx = plot( fitresult_x_L, theta_all, rRelMotion_km_L(:,1) );
% xlabel('\theta'); ylabel('x_{ref} [km]'); xlim([0,2*pi])
% set(gca,'FontSize',15);
% subplot(2,1,2)
% hz = plot( fitresult_y_L, theta_all, rRelMotion_km_L(:,2) );
% xlabel('\theta'); ylabel('y_{ref} [km]'); xlim([0,2*pi])
% set(gca,'FontSize',15);

%% 画星历相对运动分量
% figure(6)
% subplot(2,1,1)
% hx = plot( theta_DRO_eph, rv_TCR(:,1),'LineWidth',2 );
% xlabel('\theta'); ylabel('x_{eph} [km]'); xlim([0,2*pi])
% set(gca,'FontSize',15);
% subplot(2,1,2)
% hz = plot( theta_DRO_eph, rv_TCR(:,2),'LineWidth',2 );
% xlabel('\theta'); ylabel('y_{eph} [km]'); xlim([0,2*pi])
% set(gca,'FontSize',15);