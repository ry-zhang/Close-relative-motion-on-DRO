% 庞加莱截面
clear; clc;
% Poincare Section (PS) with the same engergy
% start from x axis, [x 0 0 yd]
%%
global mu
mu = 0.01215; % Earth-Moon system
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
opts_event = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events',@secYeq0);
%%
% [L] = lagrange_points(mu);
% r1 = sqrt((L(:,1)+mu).^2+L(:,2).^2);
% r2 = sqrt((L(:,1)-1+mu).^2+L(:,2).^2);
% JLp=  (L(:,1).^2+L(:,2).^2)+2*((1-mu)./r1+mu./r2); % Energy of the L-points

%%

DT = 2000; % intergration time
load('J_x0_2_3.2_.mat')
nj = 186;
aux.J = J_x0(nj,5);
aux.mu =  0.01215;
aux.dim = 6;

x0_2d = J_x0(nj,1);
% xMat = x0_2d + [-1e-5:1e-7:1e-5];
% xMat = x0_2d + [-2.3e-2:1e-4:4.4409e-2];
% xMat = x0_2d + [-2.3e-2:1e-4:0]; % 不对称
% xMat = x0_2d + [0:1e-3:4.4409e-2];
xMat = x0_2d+1e-2;
PS = cell(length(xMat), 1); % storing the PS data
Pd = cell(length(xMat), 1); % storing the periods
maxDimPd = 0;
z0 = -0.03;
% z0 = 0;
% for z0 = -0.1:1e-2:0.1
%     for z0 = -0.03
% for z0 = -0.09
isplot = 1;
% 画DRO
if isplot == 1
    figure(1)
    x0 = [x0_2d 0 0 0 ydInit2(x0_2d, 0, 0, aux) 0];
    [t0,X0] = ode113(@(t,x)odeR3bpAug(t,x,aux),[0 DT], x0, opts);
    p0 = plot3(X0(:,1),X0(:,2),X0(:,3),'LineWidth',1.5); hold on
    label('\itx \rm[LU]','\ity \rm[LU]','\itz \rm[LU]')
    box on
    axis equal
    grid on; grid minor
    xlim([0.7,1.3]); ylim([-0.5,0.5]);
end
for ii = 1:length(xMat)
    ii;
    x = xMat(ii);
    % chk - check if the initial value lies inside the zero-velocity curve
%     chk =  x^2+(2*(1-mu))/sqrt((x+mu)^2)+2*mu/sqrt((x-1+mu)^2)-J;
    temp = []; % [x y xd yd]
    T_temp = [];
    
%     if chk>0 % if good
%         yd = ydInit(x, J); % caculate the initial yd
%         xd = 0;
%         yd = ydInit1(x, xd, J); 
        X0 = [x 0 z0 0 ydInit2(x, z0, 0, aux) 0];
        
        [T,X, Te, Xe, Ie] = ode113(@(t,x)odeR3bpAug(t,x,aux),[0 DT], X0, opts_event);
        temp = [Xe Te];
        if(length(Te)>1)
            T_temp = Te(2:end)-Te(1:end-1);
            maxDimPd = max(maxDimPd, length(T_temp));
        end
%     end
    if min(X(:,1))<x0_2d-0.5
        temp = [];
    end
    PS{ii,1} = temp;
    Pd{ii,1} = T_temp;
    
    if isplot == 1 && min(X(:,1))>x0_2d-0.1
        figure(1)
        if exist('p1')
            delete(p1)
        end
        p1 = plot3(X(:,1),X(:,2),X(:,3),'Color',[217, 83, 25]/255);
    end
    
    pause(0.01)
end
hold off

pdMat = [];
for ii = 1: length(xMat)
    pdMat(:,ii) = [Pd{ii,1}; nan*ones(maxDimPd-length(Pd{ii,1}),1)];
end

%%
% close all


% hold on
for ii = 1: length(xMat)
    temp = PS{ii,1};
    if ~isempty(temp)
        figure(4);  hold on
%         plot(temp(:,1), temp(:,4),'.') ; hold on
        plot3(temp(:,1), temp(:,4), z0*ones(size(temp(:,4))),'.') ; hold on
        xlabel('\itx \rm[LU]'); ylabel('\itv_x \rm[VU]'); zlabel('\itz \rm[LU]');
        figure(5);  hold on
        plot3(temp(:,1), temp(:,3), z0*ones(size(temp(:,4))),'.') ; hold on
        xlabel('\itx \rm[LU]'); ylabel('\itz \rm[LU]'); zlabel('\itz_{\rmmax} \rm[LU]');
    end
end
% hold off
figure(4)
grid on; grid minor
xlim(x0_2d + [-0.05 0.08]);
ylim([-0.2, 0.2]);
zlim([-0.1, 0.1]);
title(['\itJ \rm= ', num2str(aux.J)])
figure(5)
grid on; grid minor
xlim(x0_2d + [-0.05 0.08]);
ylim([-0.1, 0.1]);
zlim([-0.1, 0.1]);
title(['\itJ \rm= ', num2str(aux.J)])

% figure(2)
% boxplot(pdMat,xMat)
% end