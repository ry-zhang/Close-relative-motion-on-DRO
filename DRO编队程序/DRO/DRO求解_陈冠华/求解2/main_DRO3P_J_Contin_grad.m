clear; clc
set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultTextFontSize', 14)
set(0, 'DefaultLineLineWidth', 1)

dbstop if error
%----------------------数值延拓多周期DRO轨道-------------
% 由于在延拓的收敛性不同，故分段延拓
% 使用解析梯度
% 平动点能量
%  L1--3.20033809502663
%  L2--3.184158216376
%  L3--3.02414894291943
%  L4/L5--3
%---------------------------------------------

mu = 0.01215;
aux.mu = mu;

np = 1;
x0 = 0.28;
dJ = 5e-3;
Jmat = 2.4:dJ:3.1; % continuation ok
Jmat = Jmat';

%% J=2.41延拓至J=3.0
% x0 = 0.195162398180325;
% dJ = 5e-3;
% Jmat = 2.41:dJ:3.0; % continuation ok
% Jmat = Jmat';

%%
% np = 3; %三周期轨道
aux.dim = 4; %平面问题;
xIter = x0;
dxdJ = [];

opts = odeset('RelTol',1e-13,'AbsTol',1e-13, 'Events',@secYeq0);
DT = 2*pi*10;

%%
mm = round((3.025-2.4)/dJ-1);
for ii = mm:mm+1%1:1:length(Jmat)-1
    % ii的細節并存储
    aux.J = Jmat(ii);
    yd = ydInit2(xIter, 0, aux);
    I4 = eye(4);
    X0 = [xIter, 0, 0, yd, I4(:)'];
    x0 = X0(1);  y0 = X0(2);
    vx0 = X0(3); vy0 = X0(4);
    r10 = sqrt((x0+aux.mu)^2+y0^2);
    r20= sqrt((x0-1+aux.mu)^2+y0^2);
    [~,~, Te, Xe, ~] = ode113(@(t,x)odePcr3bpAug(t,x, aux),[0 DT], X0, opts);
    XT = Xe(np+1,1:4)';  % np+1 = 回归np次
    
    Phi = Xe(np+1, 5:end);
    Phi = reshape(Phi, aux.dim, aux.dim);
    fxT = odePcr3bpAug(Te(np+1), Xe(np+1, :), aux);
    fxT = fxT(1:4);
    h = [XT(1) - X0(1); XT(3) - X0(3)];
    C = [1 0 0 0; 0 0 1 0];
    dhdx0 = -C;
    dhdxT =  C;
    dgdxT = [0 1 0 0];
    dvydx = (1/2)*(2*x0-(2*(1-aux.mu))/(x0+aux.mu)^2....
        +2*aux.mu/(1-aux.mu-x0)^2)/sqrt(x0^2....
        +(2*(1-aux.mu))/(x0+aux.mu)....
        +2*aux.mu/(1-aux.mu-x0)-aux.J);
    dvydJ =-1/(2*sqrt(x0^2+(2*(1-aux.mu))/(x0+aux.mu)+2*aux.mu/(1-aux.mu-x0)-aux.J));
    dx0dz = [1 0; 0 0; 0 0; dvydx dvydJ]; % z = [xIter J]
    dxTdx0PS = (eye(4)-fxT*dgdxT/(dgdxT*fxT))*Phi; % PS上點的變分
    dhdz = (dhdx0+dhdxT*dxTdx0PS)*dx0dz;
    
    [U,S,V] = svd(dhdz); % 延拓方向
    contDrec = V(:,2);
    dxdJ = V(1,2)/V(2,2);
    V
    if max(max(S))<1 % 分叉的處理，保持延拓方向的一致性
        dxdJOption= [V(1,1)/V(2,1), V(1,2)/V(2,2);];
        cont = abs(dxdJOption-dro3p(ii-1, 7));
        [~, id] = min(cont);
        dxdJ= V(1,id)/V(2,id);
    end
    
    dro3p(ii, :) = [x0, y0, vx0, vy0, aux.J, Te(np+1), dxdJ, S(1,1), S(2,2)];
    
    % 從ii求ii+1的初值
    aux.J = Jmat(ii+1);
    xIter = dro3p(ii,1)+dro3p(ii, 7)*dJ; % 延拓迭代初值
    yd = ydInit2(xIter, 0, aux);
    I4 = eye(4);
    X0 = [xIter, 0, 0, yd, I4(:)'];
    
    err = 1; errOld = 1;
    iterNum = 1;
    while abs(err) >1e-12%&& iterNum<1e2
        x0 = X0(1);  y0 = X0(2);
        vx0 = X0(3); vy0 = X0(4);
        r10 = sqrt((x0+aux.mu)^2+y0^2);
        r20= sqrt((x0-1+aux.mu)^2+y0^2);
        [T,X, Te, Xe, ~] = ode113(@(t,x)odePcr3bpAug(t,x, aux),[0 DT], X0, opts);
        
        XT = Xe(np+1,1:4)';  % np+1 = 回归np次
        xT = XT(1); vxT = XT(3);
        Phi = Xe(np+1, 5:end);
        Phi = reshape(Phi, aux.dim, aux.dim);
        fxT = odePcr3bpAug(Te(np+1),Xe(np+1, :), aux);
        fxT = fxT(1:4);
        
        h = [XT(1) - X0(1); XT(3) - X0(3)];
        err = norm(h);
        
        C = [1 0 0 0; 0 0 1 0];
        dhdx0 = -C;
        dhdxT =  C;
        dgdxT = [0 1 0 0];
        dvydx = (1/2)*(2*x0-(2*(1-aux.mu))/(x0+aux.mu)^2....
            +2*aux.mu/(1-aux.mu-x0)^2)/sqrt(x0^2....
            +(2*(1-aux.mu))/(x0+aux.mu)....
            +2*aux.mu/(1-aux.mu-x0)-aux.J);
        dvydJ =-1/(2*sqrt(x0^2+(2*(1-aux.mu))/(x0+aux.mu)+2*aux.mu/(1-aux.mu-x0)-aux.J));
        dx0dz = [1 0; 0 0; 0 0; dvydx dvydJ]; % z = [xIter J]
        dxTdx0PS = (eye(4)-fxT*dgdxT/(dgdxT*fxT))*Phi; % PS上點的變分
        dhdz = (dhdx0+dhdxT*dxTdx0PS)*dx0dz;
        dhdx = dhdz(:,1);
        xIter = xIter - pinv(dhdx)*h;
        
        yd = ydInit2(xIter, 0, aux);
        X0 = [xIter, 0, 0, yd, I4(:)'];
        iterNum = iterNum +1;
        if abs(errOld - err)<1e-10
            break;
        end
        if iterNum>20 % 改成最小化
            errOld = err;
        end
    end
    disp([aux.J, dro3p(ii,7),  iterNum, xIter])
end

% figure(222)
% plot(Jmat(1:length(dxdJ)),dxdJ,'-o'); hold on; grid on
% title('dxdJ')

%% 存储
% fileName = ['DRO_', num2str(np),'po_',num2str(Jmat(1)),...
%     '_',num2str(Jmat(end)),'_',num2str(dJ),'_.mat'];
% try
%     save(fileName, 'mu', 'dro3p')
% catch
%     disp('Data already exist!')
% end

%% 画图
clear;clc
load('DRO_3po_2.41_3_0.005_.mat')
mu = 0.01215;
aux.mu = mu;
aux.dim = 4;

opts = odeset('RelTol',1e-13,'AbsTol',1e-13, 'Events',@secYeq0);
close all
% figure(100)
% hold on; grid on
% axis equal
% xlabel('x [r_{EM}]')
% ylabel('y [r_{EM}]')
% data = [];
% aux.mu = mu;
% for ii = 1:10:size(dro3p, 1)
%     [T_temp, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
%     plot(X_temp(:,1),  X_temp(:,2), '-')
% end


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

%
figure(50) %% J - x圖
plot(dro3p(:,5), dro3p(:,1),'.'); hold on;  grid on
xlabel('J'); ylabel('x')
% clear