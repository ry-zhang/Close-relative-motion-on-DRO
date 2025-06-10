%% z向振动奇异点分析
% 2022-3-5
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

%% 计算t0与dt对应的phi63的值
theta1_0 = 0; theta2_0 = 0;
k0 = 0; k1 = 0; k2 = 1; 
dt = para.T0;
length_t0 = 301;
length_dt = 401;
t0_sample = linspace(0*para.T0,1*para.T0,length_t0);
dt_sample = linspace(-0.3*para.T0,3*para.T0,length_dt);
phi63_all = zeros(length_t0,length_dt);

for ii_loop = 1:length_t0
    t0 = t0_sample(ii_loop);
    rel_motion_temp_t01 = generalSol_relMotion(t0,k0,k1,k2,theta1_0,theta2_0,para,coe);
    rel_motion_temp_t02 = generalSol_relMotion(t0,k0,k1,k2,theta1_0,theta2_0+pi/2,para,coe);
    PhiSol_t0 = [rel_motion_temp_t01([3,6]),rel_motion_temp_t02([3,6])];
    
    for jj_loop = 1:length_dt
        tf = t0+dt_sample(jj_loop);
        rel_motion_temp_tf1 = generalSol_relMotion(tf,k0,k1,k2,theta1_0,theta2_0,para,coe);
        rel_motion_temp_tf2 = generalSol_relMotion(tf,k0,k1,k2,theta1_0,theta2_0+pi/2,para,coe);
        PhiSol_tf = [rel_motion_temp_tf1([3,6]),rel_motion_temp_tf2([3,6])];
        Phi = PhiSol_tf*PhiSol_t0^(-1);
        
        phi63_all(ii_loop,jj_loop) = Phi(1,2);
    end
end


%% 计算phi63的零点及其对应的phi33、phi12_tfdot的值
opts_fsolve = optimoptions('fsolve','Display','off',...
    'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'Algorithm','levenberg-marquardt');
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 para.T0], para.x0_DRO, opts);
% sol_DRO_rel = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[0 para.T0], [para.x0_DRO, zeros(1,6),reshape(eye(6,6),1,36)], opts);

[col,row] = find(phi63_all(:,1:length_dt-1)'.*phi63_all(:,2:length_dt)'<0);
phi33_all = zeros(size(row));
phi12_tfdot_all = zeros(size(row));
t0_all_zeros = t0_sample(row);
dt_all_zeros = zeros(size(row));
for kk_z0f = 1:length(col)
    t0 = t0_all_zeros(kk_z0f);
    dt_ini = dt_sample(col(kk_z0f));
    [dt,fval] = fsolve(@(x)zsingularTM(x,t0,para,coe),dt_ini,opts_fsolve);
    if abs(dt-dt_ini)>abs(2*diff(dt_sample(1:2)))
        dt_ini = dt_sample(col(kk_z0f)+1);
        [dt,fval] = fsolve(@(x)zsingularTM(x,t0,para,coe),dt_ini,opts_fsolve);
    end
    
    if max(abs(fval))>1e-5
        warning('No solution found')
    end
    dt_all_zeros(kk_z0f) = dt;
    tf = t0+dt;
    
    rel_motion_temp_t01 = generalSol_relMotion(t0,k0,k1,k2,theta1_0,theta2_0,para,coe);
    rel_motion_temp_t02 = generalSol_relMotion(t0,k0,k1,k2,theta1_0,theta2_0+pi/2,para,coe);
    PhiSol_t0 = [rel_motion_temp_t01([3,6]),rel_motion_temp_t02([3,6])];
    
    rel_motion_temp_tf1 = generalSol_relMotion(tf,k0,k1,k2,theta1_0,theta2_0,para,coe);
    rel_motion_temp_tf2 = generalSol_relMotion(tf,k0,k1,k2,theta1_0,theta2_0+pi/2,para,coe);
    PhiSol_tf = [rel_motion_temp_tf1([3,6]),rel_motion_temp_tf2([3,6])];
    Phi = PhiSol_tf*PhiSol_t0^(-1);
    phi33_all(kk_z0f) = Phi(1,1);
    
    x_DRO_t0 = deval(sol_DRO,t0)';
%     if diff([t0 tf]) ~= 0
%         sol_DRO_rel = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[t0 tf], [x_DRO_t0, zeros(1,6),reshape(eye(6,6),1,36)], opts);
%         [xphi_tf,xphidot_tf] = deval(sol_DRO_rel,tf);
%     else
%         xphidot_tf = eomM_rel3b(t0,[x_DRO_t0, zeros(1,6),reshape(eye(6,6),1,36)]',con.mu);
%     end
    I66 = eye(6,6); I66([3,6],[3,6]) = Phi;
    xphidot_tf = eomM_rel3b(t0,[x_DRO_t0, zeros(1,6),reshape(I66,1,36)]',con.mu);
    phi12_tfdot_all(kk_z0f) = xphidot_tf(45);
end

%% plot
% ---------------------------t0,dt与phi63----------------------------------
figure(1)
[t0_mesh,dt_mesh] = meshgrid(t0_sample/para.T0,dt_sample/para.T0);
surf(t0_mesh,dt_mesh,phi63_all','EdgeColor','interp')
hold on
pat = patch([min(t0_sample/para.T0),min(t0_sample/para.T0),max(t0_sample/para.T0),max(t0_sample/para.T0)],...
    [min(dt_sample/para.T0),max(dt_sample/para.T0),max(dt_sample/para.T0),min(dt_sample/para.T0)],...
    [0,0,0,0]);
set(pat,'EdgeAlpha',0,'FaceColor',0.1*[1,1,1],'FaceAlpha',0.3)
xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\Delta{\itt} [{\itT}_0]'); zlabel('$$\varphi_{12}$$','Interpreter','latex')
grid on; grid minor; box on; hold off
colorbar
xlim([min(t0_sample/para.T0),max(t0_sample/para.T0)]);
ylim([min(dt_sample/para.T0),max(dt_sample/para.T0)]);
% ylim([0,max(dt_sample/para.T0)]);
view([63,32])
% view([0,90])
% zlim([0,3])

% ---------------------------phi63的零点----------------------------------
figure(2)
[col_temp,~] = find(phi63_all(1,1:length_dt-1)'.*phi63_all(1,2:length_dt)'<0);
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
xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\Delta{\itt} [{\itT}_0]'); zlabel('$$\varphi_{36}$$','Interpreter','latex')
grid on; grid minor; box on
view([0,90])
zlim([0,3])

% ---------------------------phi63零点处的phi33----------------------------
phi33_all_reshape = reshape(phi33_all,numzeros,length(col)/numzeros);
figure(3)
colorall = get(gca,'colororder');
for ii_loop1 = 1:numzeros
    plot3(t0_all_zeros_reshape(ii_loop1,:)/para.T0,dt_all_zeros_reshape(ii_loop1,:)/para.T0,...
        phi33_all_reshape(ii_loop1,:),'LineWidth',1.5,'Color',colorall(ii_loop1,:)); hold on
    pxy = plot3(t0_all_zeros_reshape(ii_loop1,:)/para.T0,dt_all_zeros_reshape(ii_loop1,:)/para.T0,...
        0*ones(size(phi33_all_reshape(ii_loop1,:))),'LineWidth',1.5,'Color',colorall(ii_loop1,:)); hold on
    pxy.Color(4) = 0.1;
    for ii_loop2 = 1:5:length(col)/numzeros
        plot3([t0_all_zeros_reshape(ii_loop1,ii_loop2),t0_all_zeros_reshape(ii_loop1,ii_loop2)]/para.T0,...
            [dt_all_zeros_reshape(ii_loop1,ii_loop2),dt_all_zeros_reshape(ii_loop1,ii_loop2)]/para.T0,...
            [phi33_all_reshape(ii_loop1,ii_loop2),0],'LineWidth',0.5,'Color',[0.8,0.8,0.8]);
    end
end
hold off
xlim([min(t0_sample/para.T0),max(t0_sample/para.T0)]);
ylim([min(dt_sample/para.T0),max(dt_sample/para.T0)]);
% ylim([0,max(dt_sample/para.T0)]);
set(gca,'FontSize',13); 
xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\Delta{\itt} [{\itT}_0]'); zlabel('$$\varphi_{33}$$','Interpreter','latex')
grid on; grid minor; box on
set(gcf,'Renderer','painters')
view([45,13])

% ---------------------------phi63零点处的phi63导数--------------------------
phi12_tfdot_all_reshape = reshape(phi12_tfdot_all,numzeros,length(col)/numzeros);

figure(4)
colorall = get(gca,'colororder');
for ii_loop1 = 1:numzeros
    plot3(t0_all_zeros_reshape(ii_loop1,:)/para.T0,dt_all_zeros_reshape(ii_loop1,:)/para.T0,...
        phi12_tfdot_all_reshape(ii_loop1,:),'LineWidth',1.5,'Color',colorall(ii_loop1,:)); hold on
    pxy = plot3(t0_all_zeros_reshape(ii_loop1,:)/para.T0,dt_all_zeros_reshape(ii_loop1,:)/para.T0,...
        0*ones(size(phi12_tfdot_all_reshape(ii_loop1,:))),'LineWidth',1.5,'Color',colorall(ii_loop1,:)); hold on
    pxy.Color(4) = 0.1;
    for ii_loop2 = 1:5:length(col)/numzeros
        plot3([t0_all_zeros_reshape(ii_loop1,ii_loop2),t0_all_zeros_reshape(ii_loop1,ii_loop2)]/para.T0,...
            [dt_all_zeros_reshape(ii_loop1,ii_loop2),dt_all_zeros_reshape(ii_loop1,ii_loop2)]/para.T0,...
            [phi12_tfdot_all_reshape(ii_loop1,ii_loop2),0],'LineWidth',0.5,'Color',[0.8,0.8,0.8]);
    end
end
hold off
xlim([min(t0_sample/para.T0),max(t0_sample/para.T0)]);
ylim([min(dt_sample/para.T0),max(dt_sample/para.T0)]);
% ylim([0,max(dt_sample/para.T0)]);
set(gca,'FontSize',13); 
xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\Delta{\itt} [{\itT}_0]'); zlabel('$$\partial\varphi_{36}/\partial t_f$$','Interpreter','latex')
grid on; grid minor; box on
set(gcf,'Renderer','painters')
view([45,13])

%% 计算不同theta20对应的转移轨道
index_t0 = 50;
index_zeros = 5;
t0_exam = t0_all_zeros_reshape(index_zeros,index_t0);
dt_zero = dt_all_zeros_reshape(index_zeros,index_t0);
phi_exam = phi33_all_reshape(index_zeros,index_t0);
t_sample8 = linspace(0,t0_exam+dt_zero+1,2001);

% z0f = [1, 3];
z0f = [1,phi_exam];
% v0f_all = zeros(size(delta_t_all));
% theta20_all = linspace(0,2*pi,9);
theta20_all = [0,1/4,1/3,1/2,2/3,3/4]*pi;
hold off
for jj_deltat = 1:length(theta20_all)
    
    theta2_0 = theta20_all(jj_deltat);
    k2 = 1;
    theta1_0 = 0; k0 = 0; k1 = 0;
    rel_motion_temp = generalSol_relMotion(t_sample8,k0,k1,k2,theta1_0,theta2_0,para,coe);
    rel_motion_t0 = generalSol_relMotion(t0_exam,k0,k1,k2,theta1_0,theta2_0,para,coe);
    rel_motion_temp = rel_motion_temp * z0f(1)/rel_motion_t0(3,1);
%     v0f_all(jj_deltat) = sum(abs(rel_motion_temp(6,[1,end])));
    
    figure(7);
    plot(t_sample8/para.T0,rel_motion_temp(3,:),'LineWidth',1.5); hold on
    grid on; grid minor; box on
    xlim([t_sample8(1),t_sample8(end)]/para.T0);
    % ylim([-0.8,0.8]);
    set(gca,'FontSize',13)
    xlabel('{\itt} [{\itT}_0]'); ylabel('\itz_L');
end

xlim([t_sample8(1),t_sample8(end)]/para.T0);
str = strings(1,length(theta20_all));
for jj_theta20 = 1:length(theta20_all)
    delta_t = theta20_all(jj_theta20);
    str(jj_theta20) = ['\theta_{2,0}=',num2str(delta_t/pi,'%.2f'),'\pi'];
end
% legend(str,'Location','north')
legend(str,'Location','northeastoutside','AutoUpdate', 'off')
xticks([0,1,2,3])

%% 计算phi12对tf的导数
% index_t0 = 10;
% index_zeros = 2;
% t0_exam = t0_all_zeros_reshape(index_zeros,index_t0);
% dt_exam = dt_all_zeros_reshape(index_zeros,index_t0);
% phi_exam = phi33_all_reshape(index_zeros,index_t0);
% hold off
% theta2_0_all = linspace(0,2*pi,201);
% ratio_all = zeros(size(theta2_0_all));
% isplot = 0;
% for jj_theta2 = 1:length(theta2_0_all)
%     theta2_0 = theta2_0_all(jj_theta2);
%     
%     t_sample9 = linspace(t0_exam,t0_exam+dt_exam,100);
%     theta1_0 = 0; k0 = 0; k1 = 0; k2 = 1;
%     rel_motion_temp = generalSol_relMotion(t_sample9,k0,k1,k2,theta1_0,theta2_0,para,coe);
%     ratio_all(jj_theta2) = rel_motion_temp(3,end)/rel_motion_temp(3,1);
%     if isplot == 1
%         figure(8);
%         plot(t_sample9/para.T0,rel_motion_temp(3,:),'LineWidth',1.5); hold on
%         grid on; grid minor; box on
%         xlim([t_sample9(1),t_sample9(end)]/para.T0);
%         % ylim([-0.8,0.8]);
%         set(gca,'FontSize',13)
%         xlabel('{\itt}_0 [{\itT}_0]'); ylabel('\itz_L');
%     end
% end
% figure(8)
% plot(theta2_0_all,ratio_all,'LineWidth',1.5);
% grid on; grid minor; box on
% xlim([0,2*pi]);
% ylim(max(ratio_all)+[-0.01,0.01]);
% set(gca,'FontSize',13)
% xlabel('\theta_2(0) [rad]'); ylabel('z(t_f)/z(t_0)');
