% 延拓求解
clear; clc
set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultTextFontSize', 14)
set(0, 'DefaultLineLineWidth', 2)

dbstop if error
%----------------------数值延拓多周期DRO轨道-------------
% 由于在延拓的收敛性不同，故分段延拓
% 平动点能量
%  L1--3.20033809502663
%  L2--3.184158216376
%  L3--3.02414894291943
%  L4/L5--3
%---------------------------------------------

mu = 0.01215;
aux.mu = mu;

np = 1; 

%% 1
x0 = 0.123868449285141;
dJ = 1e-3;
Jmat = 2:dJ:3.2;

%% 3
% x0 = 0.082302039412579;
% dJ = -5e-3;
% Jmat = 2.125:dJ:2;

%% 6周期
% x0 = 0.93739;
% x0 = 0.885651463314672;
% dJ = -2e-3;
% Jmat = 3.06:dJ:3.05

%% 5周期
% % x0 = 0.935236889218813;
% x0 = 0.935095754388488;
% % x0 = 0.870760142884322;
% 
% x0 = 0.2315;
% x0 = 0.326;
% x0 = 0.95465;
% dJ = 2e-3;
% Jmat = 3.06:dJ:3.1


%% 第一段 J=2.7逆向延拓至J=2.4
% x0 = 0.416884406153494;
% dJ = -5e-3;
% Jmat = 2.7:dJ:2.4;  % continuation ok
% % Jmat = Jmat';

%% 第二段 J=2.7正向延拓至J=3.0
% x0 = 0.416884406153494;
% dJ = 5e-3;
% Jmat = 2.7:dJ:3;

%% 第三段 J=3.0正向延拓至J=3.1
% x0 = 0.831641070677090;
% dJ = 2e-3;
% Jmat = 3:dJ:3.1;
% Jmat = Jmat';

%%% 第二段 J=3.1正向延拓至J=3.142
% x0 = [0.820597513226115];
% dJ = 2e-3;
% Jmat = 3.1:dJ:3.142;
% Jmat = Jmat';
% 
% 
% x0 = [0.819759333449791];
% dJ = 1e-3;
% Jmat = 3.142:dJ:3.154;
% Jmat = Jmat';

%%
% np = 3; %三周期轨道

xIter = x0;
dxdJ = [(0.082302039412579-0.083775507907209)/dJ];
for ii = 2:1:length(Jmat)
    ii
    J = Jmat(ii);
    aux.J = J;
    yd = ydInit1(xIter, 0, aux);
    X0 = [xIter, 0, 0, yd];
    
    DT = 2*pi*10;   
%     DT = 5*pi*10;  
    opts = odeset('RelTol',1e-7,'AbsTol',1e-7, 'Events',@secYeq0);
    
    err = 1;
    iterNum = 1;
    if ~isempty(dxdJ)
        xIter = xIter +dJ*dxdJ(end);
    else
         xIter = xIter +dJ;
    end
    
    while abs(err) >1e-12 %&& iterNum<1e2
        [~,~, Te, Xe, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 DT], X0, opts);        
        xT = Xe(np+1,1); xdT = Xe(np+1,3);       % 4 = 回归三次
        dx = 1e-6;        
        X0Var = [xIter+dx, 0, 0, ydInit1(xIter+dx, 0, aux)];
        [~,~, TeVar, XeVar, ~] = ode113(@(t,x)odePcr3bp(t ,x, aux),[0 DT], X0Var, opts);
        try
            p_xTx0 = (XeVar(np+1,1)-Xe(np+1,1))/dx;
        catch
            [~,~, TeVar, XeVar, ~] = ode113(@(t,x)odePcr3bp(t ,x, aux),[0 2*DT], X0Var, opts);
            p_xTx0 = (XeVar(np+1,1)-Xe(np+1,1))/dx;
        end
        p_xdTx0 = (XeVar(np+1,3)-Xe(np+1,3))/dx;
        
        if 1
            err = ((xT-xIter)^2+xdT^2)/2;
            Jac = (xT-xIter)*(p_xTx0-1)+xdT*p_xdTx0;
        else
            err = xT-xIter;
            Jac = p_xTx0-1;
        end        
        
        xIter = xIter - err/Jac;
        yd = ydInit1(xIter, 0, aux);
        X0 = [xIter, 0, 0, yd];
        iterNum = iterNum +1;
    end
    
    x0(ii,1) = xIter;
    dxdJ(ii-1,1) = (x0(ii,1)- x0(ii-1,1))/dJ;
    [J, iterNum, xIter, dxdJ(ii-1,1)] 
end

figure(101)
plot(Jmat(1:length(dxdJ)),dxdJ,'-o'); hold on; grid on
title('dxdJ')
%% 存储
for ii = 1:length(x0)    
    aux.J = Jmat(ii);
    yd = ydInit1(x0(ii), 0, aux);
    
    [~, ~, Te, Xe, ~] = ode113(@(t,x)odePcr3bp(t ,x, aux),[0 DT], [x0(ii), 0, 0, yd ], opts);   
    % posvel, J ,T
    dro3p(ii, :) = [x0(ii), 0, 0, yd, Jmat(ii), Te(np+1)];
end

fileName = ['DRO_', num2str(np),'po_',num2str(Jmat(1)),...
    '_',num2str(Jmat(end)),'_',num2str(dJ),'_.mat'];
try
    save(fileName, 'mu', 'dro3p')
catch
    disp('Data already exist!')
end



%% 画图
opts = odeset('RelTol',1e-10,'AbsTol',1e-10, 'Events',@secYeq0);
close all
figure(100)
hold on; grid on
axis equal
xlabel('x [r_{EM}]')
ylabel('y [r_{EM}]')
data = [];

for ii = 1:size(dro3p, 1)
    [T_temp, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
    plot(X_temp(:,1), X_temp(:,2),'-')
end
% for ii = size(dro3p, 1)-1  %2:size(dro3p, 1)    
%     [T_temp, X_temp, ~, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 dro3p(ii, 6)],  dro3p(ii, 1:4), opts);
%     plot(X_temp(:,1), X_temp(:,2),'-')
% end

