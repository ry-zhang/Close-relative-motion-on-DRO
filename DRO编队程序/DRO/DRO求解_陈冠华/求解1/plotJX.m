% 画画专用
% clear all; clc

%% 画所有
% load('J_x0_2_3.2_.mat')
% load('J_x3_2.125_3.1_.mat')
% load('J_x4a_2_2.4_2.75_.mat')
% load('J_x4b_2_2.4_2.75_.mat')
% load('J_x4c_2.995_3.1_.mat')
% load('J_x4d_2.994_3.06_.mat')
% load('J_x5a_2.4_2.45_.mat')
% load('J_x5b_2.4_2.45_.mat')
% load('J_x5c_3.026_3.1_.mat')
% load('J_x5d_3.025_3.09_.mat')
% load('J_x6a_3.05_3.1_.mat')
% load('J_x6b_3.05_3.1_.mat')
% 
% 
% % figure(1)
% % plot(J_x3(:,5),J_x3(:,1),...
% %     J_x4_a(:,5),J_x4_a(:,1),...
% %     J_x4_b(:,5),J_x4_b(:,1),...
% %     J_x4_c(:,5),J_x4_c(:,1),...     
% %     J_x4_d(:,5),J_x4_d(:,1),...
% %     J_x5_a(:,5),J_x5_a(:,1),...
% %     J_x5_b(:,5),J_x5_b(:,1),...
% %     J_x5_c(:,5),J_x5_c(:,1),... 
% %     J_x5_d(:,5),J_x5_d(:,1),...
% %     J_x6_a(:,5),J_x6_a(:,1),...
% %     J_x6_b(:,5),J_x6_b(:,1),...
% %     J_x0(:,5),J_x0(:,1),...
% %     'LineWidth',1.5);
% % grid on;
% figure(1)
% f3 = plot(J_x3(:,5),J_x3(:,1),'LineWidth',1.5);hold on
% f4 = plot(J_x4_a(:,5),J_x4_a(:,1),'LineWidth',1.5);hold on
% plot(J_x4_b(:,5),J_x4_b(:,1),...
%     J_x4_c(:,5),J_x4_c(:,1),...     
%     J_x0(:,5),J_x0(:,1),'LineWidth',1.5);hold on
% f5 = plot(J_x5_a(:,5),J_x5_a(:,1),'LineWidth',1.5);hold on
% plot(J_x5_b(:,5),J_x5_b(:,1),...
%     J_x5_c(:,5),J_x5_c(:,1),... 
%     J_x5_d(:,5),J_x5_d(:,1),'LineWidth',1.5);hold on
% f6 = plot(J_x6_a(:,5),J_x6_a(:,1),'LineWidth',1.5);hold on
% plot(J_x6_b(:,5),J_x6_b(:,1),'LineWidth',1.5);hold on
% f1 = plot(J_x0(:,5),J_x0(:,1),'LineWidth',1.5);hold on
% 
% xlabel('J');ylabel('x0')
% xlim([2.2,3.1]);
% % legend('3-period','4-period','4-period','4-period',...
% %     '4-period','5-period','5-period','5-period','5-period',...
% %     '6-period','6-period','1-period')
% legend([f3,f4,f5,f6,f1],'3-period','4-period','5-period','6-period','1-period')
% hold on;

%% 3d
clear all;clc
% load('JXT_3d_.mat')
% load('JXT_3d_5pJ=2.99_.mat')
% load('J_x5c_3.026_3.1_.mat')
% load('J_x0_2_3.2_0.001_.mat')
% load('JXT_3d_5pJ=2.99_.mat')
% load('JXT_3d_5pJ=2.9_.mat')
JX_3d = JX_3d_5p;

global mu
mu = 0.01215;
aux.mu = mu;
aux.dim = 6;


figure
hold on;
% axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%  1738/384401=0.0045;
rm = 0.0045;
[xm,ym,zm] = sphere(40);
Xm = xm*rm+1;
Ym = ym*rm;
Zm = zm*rm;
% surf(Xm,Ym,Zm);
theta=0:pi/100:2*pi;
xm_ = 1+rm*cos(theta);
ym_ = rm*sin(theta);
zm_ = rm*cos(theta);

% plot3(0,0,0,'square');hold on
% plot3(1,0,0,'square');hold on

opts = odeset('RelTol',1e-13,'AbsTol',1e-13, 'Events',@secYeq0);

% z0 = 0.02*ones(length(J_x5_c),1);
% JX_3d = [J_x5_c(:,1:2), z0, J_x5_c(:,3:4), zeros(length(J_x5_c),1), J_x5_c(:,5:6)];
% JX_3d(:,3) = z0;

figure(300)
for ii = 10:20:length(JX_3d(:,1))
    
    DT = 60;
[T_temp, X_temp, ~, ~] = ode113(@(t,x)odeR3bpAug(t,x, aux),[0 DT],  JX_3d(ii,1:6), opts);
X = [ X_temp(:,1), X_temp(:,2), X_temp(:,3), JX_3d(ii,7)*ones(length(X_temp(:,1)),1)];
% plot3(X(:,1),X(:,2),X(:,3),'-','linewidth',1.2);hold on
% X(end,3) = NaN;
% patch(X(:,1),X(:,2),X(:,3),X(:,4),'EdgeColor','interp','MarkerFaceColor','flat');hold on

%
subplot(2,2,1); 
% plot3(0,0,0,'square');hold on
plot3(X(:,1),X(:,2),X(:,3),'-','linewidth',1.2,'color',[.7 .7 .7]);hold on; grid on; axis equal
surf(Xm,Ym,Zm);
xlabel('x'); ylabel('y')
subplot(2,2,2); 
plot(X(:,1),X(:,2), '-','linewidth',1.2,'color',[.7 .7 .7]); hold on; grid on; axis equal
plot(xm_,ym_);hold on
xlabel('x'); ylabel('y')
subplot(2,2,3); 
plot(X(:,1),X(:,3), '-','linewidth',1.2,'color',[.7 .7 .7]); hold on; grid on; axis equal
plot(xm_,ym_);hold on
xlabel('x'); ylabel('z')
subplot(2,2,4); 
plot(zm_,ym_);hold on
plot(X(:,2),X(:,3), '-','linewidth',1.2,'color',[.7 .7 .7]); hold on; grid on; axis equal
xlabel('y'); ylabel('z')

end

% colorbar
title(['J≥2.99'])
% ylabel(colorbar,'Jacobi energy');
% view([150,20])

% plot(X(:,1), X(:,2),'-','linewidth',1.5)

%% 2d
clear all;clc
% load('J_x0_2_3.2_0.001_.mat')
% load('J_x4a_2_2.4_2.72_.mat')
% load('J_x4c_2.996_3.076_.mat')
% load('J_x3_2_2.4_3_.mat')
load('J_x5a_2.4_2.45_.mat')
J_x0 = J_x5_a;

global mu
mu = 0.01215;
aux.mu = mu;
aux.dim = 4;

figure
hold on;
% axis equal

% plot(0,0,'square');hold on
plot(1,0,'square');hold on

opts = odeset('RelTol',1e-13,'AbsTol',1e-13, 'Events',@secYeq0);

for ii = 1:length(J_x0(:,1))-5
    
%     DT = 20;
[T_temp, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 J_x0(ii,6)],  J_x0(ii,1:4), opts);
X = [ X_temp(:,1), X_temp(:,2),  J_x0(ii,5)*ones(length(X_temp(:,1)),1)];
% plot(X(:,1), X(:,2),'-');hold on
% plot(X(:,1), X(:,2),'-','linewidth',1.2);hold on
X(end,3) = NaN;
patch(X(:,1),X(:,2),X(:,3),'EdgeColor','interp','MarkerFaceColor','flat');hold on
% plot3(X(:,1),X(:,2),X(:,3),'-','linewidth',1.2,'color',[.7 .7 .7]);hold on

end
plot(X(:,1), X(:,2),'-','linewidth',1.2);hold on
axis equal
colorbar
xlabel('x','FontName','Times New Roman','FontSize',17)
ylabel('y','FontName','Times New Roman','FontSize',17)
zlabel('z','FontName','Times New Roman','FontSize',17)
ylabel(colorbar,'Jacobi energy','FontName','Times New Roman','FontSize',17);
% colorbar;

% plot(X(:,1), X(:,2),'-','linewidth',1.5)

%% subplot
clear;clc
load('DRO_3po_2.41_3_0.005_.mat')
mu = 0.01215;
aux.mu = mu;
aux.dim = 4;
opts = odeset('RelTol',1e-13,'AbsTol',1e-13, 'Events',@secYeq0);
close all

figure(101) %% 分叉檢查
plot(dro3p(:,5), dro3p(:,8),'-o'); hold on; grid on
plot(dro3p(:,5), dro3p(:,9),'-o'); hold on; grid on


figure(102) %% 一系列軌道圖
subplot(2,3,1); ii = 1;
[~, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
plot(X_temp(:,1),  X_temp(:,2), '-'); hold on; grid on; axis equal
xlabel('x'); ylabel('y')
title(['J=',num2str(dro3p(ii, 5))]);

subplot(2,3,2); ii = 89; axis equal
[~, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
plot(X_temp(:,1),  X_temp(:,2), '-'); hold on; grid on; axis equal
xlabel('x'); ylabel('y')
title(['J=',num2str(dro3p(ii, 5))]);

subplot(2,3,3); ii = 99; axis equal
[~, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
plot(X_temp(:,1),  X_temp(:,2), '-'); hold on; grid on; axis equal
xlabel('x'); ylabel('y')
title(['J=',num2str(dro3p(ii, 5))]);

subplot(2,3,4); ii = 110; axis equal
[~, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
plot(X_temp(:,1),  X_temp(:,2), '-'); hold on; grid on; axis equal
xlabel('x'); ylabel('y')
title(['J=',num2str(dro3p(ii, 5))]);

subplot(2,3,5); ii = 118; axis equal
[~, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
plot(X_temp(:,1),  X_temp(:,2), '-'); hold on; grid on; axis equal
xlabel('x'); ylabel('y')
title(['J=',num2str(dro3p(ii, 5))]);

%% J-T
clear all;clc
load('J_x0_2_3.2_0.001_.mat')
load('J_x3_2.4_3.1_.mat')
load('J_x4a_2_2.4_2.75_.mat')
load('J_x4c_2.995_3.1_.mat')
load('J_x5a_2.4_2.45_.mat')
load('J_x5c_3.026_3.1_.mat')
load('JXT_3d_.mat')
load('JXT_3d_5pJ=2.99_.mat')
load('JXT_3d_5pJ=2.9_.mat')

figure
plot(J_x3(:,end-1),J_x3(:,end)/3,...
    J_x4_a(:,end-1),J_x4_a(:,end)/4,...
    J_x5_a(:,end-1),J_x5_a(:,end)/5,...
    JX_3d(:,end-1),JX_3d(:,end),...
    J_x4_c(:,end-1),J_x4_c(:,end)/4,...
    J_x5_c(:,end-1),J_x5_c(:,end)/5,...
    'LineWidth',1.5);
xlim([2.4,3.1])
legend('P3','P4','P5','DRO')
xlabel('J');ylabel('T/N')

%%

 load('JXT_3d_5pJ=2.99_.mat')
dro3p = JX_3d_5p;

figure(50) %% J - x圖
plot(dro3p(:,8), dro3p(:,1),'.'); hold on;  grid on
xlabel('J'); ylabel('x')
title('3d-P5DRO')

load('data_3dP5DRO_.mat')
dro3p = data3dP5DRO;
figure(50) %% J - x圖
plot(dro3p(:,8), dro3p(:,3),'.'); hold on;  grid on
xlabel('J'); ylabel('z')
title('3d-P5DRO')

% 3d
load('JXT_3d_.mat')
dro3p = JX_3d;

figure(50) %% J - x圖
plot(dro3p(:,7), dro3p(:,3),'-'); hold on;  grid on
xlabel('J'); ylabel('|z|')
title('3d-DRO')
