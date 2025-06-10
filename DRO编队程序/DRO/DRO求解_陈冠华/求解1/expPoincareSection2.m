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

DT = 1000; % intergration time
load('J_x0_2_3.2_.mat')
nj = 189;
aux.J = J_x0(nj,5);
aux.mu =  0.01215;
aux.dim = 6;
x0_2d = J_x0(nj,1);

aux.J = 2.93052374520898;
x0_2d = 0.808936204314423;
% aux.J = 2.93584;
% x0_2d = 0.81255;

% xMat = x0_2d + [-1e-5:1e-7:1e-5];
% xMat = x0_2d + [-3e-2:1e-3:3e-2];
% xMat = x0_2d + [-2.3e-2:1e-4:0]; % 不对称
xMat = x0_2d + [0:1e-4:4.4409e-2];
% xMat = x0_2d + [0.00885:1e-4:0.02895];
% xMat = [0.831636204314423:1e-5:0.831936204314423]; %11:4
% xMat = x0_2d + [0:1e-4:1.4409e-2]+1e-2;
% xMat = x0_2d + [0.021:1e-3:0.038];
% xMat = x0_2d+1e-2;

z0_all = -0.09:1e-2:0.09;
% z0_all = -0.03;

PS = cell(length(z0_all),length(xMat), 1); % storing the PS data
Pd = cell(length(z0_all),length(xMat), 1); % storing the periods
maxDimPd = 0;
% for jj = 1:length(z0_all)
for jj = 10
    z0 = z0_all(jj);
    isplot = 0;
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
%         view(0,0)
        view(0,90)
        xlim([0.7,1.3]); ylim([-0.5,0.5]);
    end

%     for ii = 1:length(xMat)
%     for ii = 24
    parfor ii = 1:length(xMat)
        ii;
        x = xMat(ii);
%         x = 0.831806;
%         x = 0.840601204314423;
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
%             sol = ode113(@(t,x)odeR3bpAug(t,x,aux),[0 DT], X0, opts);
%             t_sample = linspace(0,DT,1e5);
%             X_sample = deval(sol,t_sample);
            temp = [Xe Te];
            if(length(Te)>1)
                T_temp = Te(2:end)-Te(1:end-1);
                maxDimPd = max(maxDimPd, length(T_temp));
            end
    %     end
        if min(X(:,1))<x0_2d-0.1
            temp = [];
        end
        PS{jj,ii,1} = temp;
        Pd{jj,ii,1} = T_temp;

        if isplot == 1 && min(X(:,1))>x0_2d-0.1
            figure(1)
            if exist('p1')
                delete(p1)
            end
            p1 = plot3(X(:,1),X(:,2),X(:,3),'Color',[217, 83, 25]/255);
            p1.Color(4) = 0.5;
        end

        pause(0.01)
    end
    hold off

    pdMat = [];
    for ii = 1: length(xMat)
        pdMat(:,ii) = [Pd{jj,ii,1}; nan*ones(maxDimPd-length(Pd{jj,ii,1}),1)];
    end

end

    
%% 画庞加莱截面
% close all
for jj = 1:length(z0_all)
    z0 = z0_all(jj);
    for ii = 1: length(xMat)
        temp = PS{jj,ii,1};
        if ~isempty(temp)
            figure(2);  hold on
    %         plot(temp(:,1), temp(:,4),'.') ; hold on
            plot3(temp(:,1), temp(:,4), z0*ones(size(temp(:,4))),'.') ; hold on
%             figure(3);  hold on
%             plot3(temp(:,1), temp(:,3), z0*ones(size(temp(:,4))),'.') ; hold on
        end
    end
end
% hold off
figure(2)
grid on; grid minor
xlim(x0_2d + [-0.05 0.08]);
ylim([-0.2, 0.2]);
zlim([-0.1, 0.1]);
title(['\itJ \rm= ', num2str(aux.J)])
xlabel('\itx \rm[LU]'); ylabel('\itv_x \rm[VU]'); zlabel('\itz_{\rmmax} \rm[LU]');

% figure(3)
% grid on; grid minor
% xlim(x0_2d + [-0.05 0.08]);
% ylim([-0.1, 0.1]);
% zlim([-0.1, 0.1]);
% title(['\itJ \rm= ', num2str(aux.J)])
% xlabel('\itx \rm[LU]'); ylabel('\itz \rm[LU]'); zlabel('\itz_{\rmmax} \rm[LU]');

% save a_expPonSec aux z0_all xMat PS Pd