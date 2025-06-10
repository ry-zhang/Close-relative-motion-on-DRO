%% xy方向振动奇异点分析
% 2022-3-9
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
clear

format longg
format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 13);%坐标轴
set(0,'defaultTextFontSize', 13);%文字
set(0,'defaultLineLineWidth',1.5)

addpath('../../subF_eom(CR3BP)')
load('generalSolFFT_12.mat')

%% 计算t0与dt对应的phi21的值
length_t0 = 101;
length_dt = 1001;
t0_sample = linspace(0*para.T0,1*para.T0,length_t0);
dt_sample = linspace(0.01,3*para.T0,length_dt);
phidet_all = zeros(length_t0,length_dt);

opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 para.T0], para.x0_DRO, opts);

parfor ii_loop = 1:length_t0
    t0 = t0_sample(ii_loop);
    xt0_DRO = deval(sol_DRO,t0)';
    sol_DRO_rel = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[t0 4*para.T0], [xt0_DRO, zeros(1,6),reshape(eye(6,6),1,36)], opts);
    for jj_loop = 1:length_dt
        tf = t0+dt_sample(jj_loop);
        [xphi_tf,xphidot_tf] = deval(sol_DRO_rel,tf);
        Phi = reshape(xphi_tf(13:end),6,6);
        phidet_all(ii_loop,jj_loop) = det(Phi([1,2],[4,5]));
    end
end

%% 计算det(Phi([1,2],[4,5]))的零点及其对应的相位值
opts_fsolve = optimoptions('fsolve','Display','off',...
    'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'Algorithm','levenberg-marquardt');
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
[col,row] = find(phidet_all(:,1:length_dt-1)'.*phidet_all(:,2:length_dt)'<0);
t0_all_zeros = t0_sample(row);
dt_all_zeros = zeros(size(row));
k_all = zeros(size(row));
% phi12_tfdot_all = zeros(size(row));
parfor kk_z0f = 1:length(col)
    t0 = t0_all_zeros(kk_z0f);
    dt_ini = dt_sample(col(kk_z0f));
    xt0_DRO = deval(sol_DRO,t0);
    [dt,fval] = fsolve(@(x)xysingularTM(x,t0,xt0_DRO,con),dt_ini,opts_fsolve);
    if abs(dt-dt_ini)>abs(2*diff(dt_sample(1:2)))
        dt_ini = dt_sample(col(kk_z0f)+1);
        [dt,fval] = fsolve(@(x)xysingularTM(x,t0,para,coe),dt_ini,opts_fsolve);
    end
    
    tf = t0+dt;
    [~,sol_DRO_rel] = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[t0 tf], [xt0_DRO', zeros(1,6),reshape(eye(6,6),1,36)], opts);
    Phi = reshape(sol_DRO_rel(end,13:end),6,6);
    Phiv = (Phi(1:2,4:5));
    k_all(kk_z0f) = Phiv(2,2)/Phiv(1,2);
    
    if max(abs(fval))>1e-5
        warning('No solution found')
    end
    
    dt_all_zeros(kk_z0f) = dt;
    
%     tf = t0+dt;
%     x_DRO_t0 = deval(sol_DRO,t0)';
%     if diff([t0 tf]) ~= 0
%         sol_DRO_rel = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[t0 tf], [x_DRO_t0, zeros(1,6),reshape(eye(6,6),1,36)], opts);
%         [xphi_tf,xphidot_tf] = deval(sol_DRO_rel,tf);
%     else
%         xphidot_tf = eomM_rel3b(t0,[x_DRO_t0, zeros(1,6),reshape(eye(6,6),1,36)]',con.mu);
%     end
%     I66 = eye(6,6); I66([3,6],[3,6]) = Phi;
%     xphidot_tf = eomM_rel3b(t0,[x_DRO_t0, zeros(1,6),reshape(I66,1,36)]',con.mu);
%     phi12_tfdot_all(kk_z0f) = xphidot_tf(45);
end

%% plot
% ---------------------------t0,dt与||\Phi_{xy}^{v}||----------------------------------
figure(1)
[t0_mesh,dt_mesh] = meshgrid(t0_sample/para.T0,dt_sample/para.T0);
surf(t0_mesh,dt_mesh,phidet_all','EdgeColor','interp')
hold on
pat = patch([min(t0_sample/para.T0),min(t0_sample/para.T0),max(t0_sample/para.T0),max(t0_sample/para.T0)],...
    [min(dt_sample/para.T0),max(dt_sample/para.T0),max(dt_sample/para.T0),min(dt_sample/para.T0)],...
    [0,0,0,0]);
set(pat,'EdgeAlpha',0,'FaceColor',0.1*[1,1,1],'FaceAlpha',0.3)
xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\Delta{\itt} [{\itT}_0]'); zlabel('$$||\Phi_{xy}^{v}||$$','Interpreter','latex')
grid on; grid minor; box on; hold off
colorbar
xlim([min(t0_sample/para.T0),max(t0_sample/para.T0)]);
ylim([min(dt_sample/para.T0),max(dt_sample/para.T0)]);
% ylim([0,max(dt_sample/para.T0)]);
view([72,21])

% ---------------------------||\Phi_{xy}^{v}||的零点----------------------------------
figure(2)
[col_temp,~] = find(phidet_all(1,1:length_dt-1)'.*phidet_all(1,2:length_dt)'<0);
numzeros = length(col_temp);
t0_all_zeros_reshape = reshape(t0_all_zeros,numzeros,length(col)/numzeros);
dt_all_zeros_reshape = reshape(dt_all_zeros,numzeros,length(col)/numzeros);
for ii_loop1 = 1:numzeros
    plot(t0_all_zeros_reshape(ii_loop1,:)/para.T0,dt_all_zeros_reshape(ii_loop1,:)/para.T0,...
        'LineWidth',1.5,'Color',[0, 114, 189]/255); hold on
end
ylim([min(dt_sample/para.T0),max(dt_sample/para.T0)]);
% plot([1,0.5],[0,1])
% plot([0,0.5],[1,0]) 
hold off
xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\Delta{\itt} [{\itT}_0]');
grid on; grid minor; box on
view([0,90])

% ---------------------------||\Phi_{xy}^{v}||零点处的k----------------------------
figure(3)
k_all_reshape = reshape(k_all,numzeros,length(col)/numzeros);
colorall = get(gca,'colororder');
for ii_loop1 = 1:numzeros
    plot3(t0_all_zeros_reshape(ii_loop1,:)/para.T0,dt_all_zeros_reshape(ii_loop1,:)/para.T0,...
        k_all_reshape(ii_loop1,:),'.','LineWidth',1.5,'Color',colorall(ii_loop1,:)); hold on
    pxy = plot3(t0_all_zeros_reshape(ii_loop1,:)/para.T0,dt_all_zeros_reshape(ii_loop1,:)/para.T0,...
        0*ones(size(k_all_reshape(ii_loop1,:))),'LineWidth',1.5,'Color',colorall(ii_loop1,:)); hold on
    pxy.Color(4) = 0.1;
    for ii_loop2 = 1:5:length(col)/numzeros
        plot3([t0_all_zeros_reshape(ii_loop1,ii_loop2),t0_all_zeros_reshape(ii_loop1,ii_loop2)]/para.T0,...
            [dt_all_zeros_reshape(ii_loop1,ii_loop2),dt_all_zeros_reshape(ii_loop1,ii_loop2)]/para.T0,...
            [k_all_reshape(ii_loop1,ii_loop2),0],'LineWidth',0.5,'Color',[0.8,0.8,0.8]);
    end
end
hold off
xlim([min(t0_sample/para.T0),max(t0_sample/para.T0)]);
ylim([min(dt_sample/para.T0),max(dt_sample/para.T0)]);
% ylim([0,max(dt_sample/para.T0)]);
set(gca,'FontSize',13); 
xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\Delta{\itt} [{\itT}_0]'); zlabel('\itk_{xy}')
grid on; grid minor; box on
set(gcf,'Renderer','painters')
view([45,13])

%% 积分计算燃耗
t0 = 0.1*para.T0;
xt0_DRO = deval(sol_DRO,t0);

x0_REL = [-5,5,0,0,0,0]'/con.r_norma;
rf_rel = [10,10,0]'/con.r_norma;
x0 = [0,0,0]';
dv_all = zeros(length_dt,3);
fval_all = zeros(length_dt,3);
parfor jj_loop = 1:length_dt
    dt = dt_sample(jj_loop);
    
    % 线性化模型下，用状态转移矩阵求解变轨，因此线性化模型下的变轨问题有且仅有唯一解
    [~,sol_DRO_rel] = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[t0 t0+dt], [xt0_DRO', zeros(1,6),reshape(eye(6,6),1,36)], opts);
    Phi = reshape(sol_DRO_rel(end,13:end),6,6);
    Phiv = Phi(1:3,4:6);
    Phir = Phi(1:3,1:3);
    v0 = Phiv^(-1)*(rf_rel - Phir*x0_REL(1:3));
    dv_all(jj_loop,:) = v0';
    % 线性化模型下，用fsolve迭代求解变轨
%     [x,fval] = fsolve(@(x)maneuver_relmotion(x,dt,xt0_DRO,x0_REL,rf_rel,con),x0,opts_fsolve);
%     x_all(jj_loop,:) = x';
%     fval_all(jj_loop,:) = fval';
    
    % 二者误差
%     error_all(jj_loop,:) = v0'-x';
end
dv_norm_all = sqrt(sum(dv_all.^2,2))*con.v_norma*1e3;

% ---------------------------计算变轨脉冲----------------------------------
figure(7)
semilogy(dt_sample/para.T0,dv_norm_all)
xlabel('\Delta{\itt} [{\itT}_0]'); ylabel('\Delta{\itv} [m/s]');
grid on; grid minor; box on
view([0,90])
% ---------------------------画奇异点变轨轨道----------------------------------
[~,index_all] = findpeaks(dv_norm_all);
index = index_all(2)-200;
[~,sol_DRO_rel] = ode113(@(t,x)eom_rel3b(t,x,con.mu),[t0 t0+dt_sample(index)], [xt0_DRO; x0_REL+[zeros(1,3),dv_all(index,:)]'], opts);
figure(8)
plot3(sol_DRO_rel(:,7)*con.r_norma,sol_DRO_rel(:,8)*con.r_norma,sol_DRO_rel(:,9)*con.r_norma)
axis equal
xlabel('{\itx} [km]'); ylabel('{\ity} [km]'); zlabel('{\itz} [km]');
grid on; grid minor; box on
view([0,90])