% function expFind_Dro_J_Contu
clear; clc

mu = 0.01215;
aux.mu = mu;

x0 = 0.284650774667595;
dJ = 1e-1;
Jmat = 2.4:dJ:3.4; % continuation ok
Jmat = Jmat';

xIter = x0;
dxdJ = [];
for ii = 6:1:length(Jmat)
    J = Jmat(ii);
    aux.J = J;
    yd = ydInit2(xIter, 0, aux);
    X0Iter = [xIter, 0, 0, yd];
    
    DT = 2*pi*10;
    opts = odeset('RelTol',1e-13,'AbsTol',1e-13, 'Events',@secYeq0);
    
    err = 1;
    iterNum = 1;
    if ~isempty(dxdJ)
        xIter = xIter +dJ*dxdJ(end);
    else
        xIter = xIter +dJ;
    end
    
    while abs(err) >1e-12 %&& iterNum<1e2
        [t_temp,X_temp, Te, Xe, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 DT], X0Iter, opts);
        
        xT = Xe(4,1); xdT = Xe(4,3); 
        
        dx = 1e-6; % used to caculate  p_xTx0
        
        X0Var = [xIter+dx, 0, 0, ydInit2(xIter+dx, 0, aux)];
        
        [~,~, ~, XeVar, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 DT], X0Var, opts);
        try
            p_XTx0 = (XeVar(4,:)-Xe(4,:))/dx;
        catch
            [~,~, ~, XeVar, ~] = ode113(@(t,x)odePcr3bp(t,x, aux),[0 2*DT], X0Var, opts);
            p_XTx0 = (XeVar(4,:)-Xe(4,:))/dx;
        end
        
        p_xTx0 = p_XTx0(1);
        p_xdTx0 = p_XTx0(3);
        
        err = ((xT-xIter)^2+xdT^2)/2;
        Jac = (xT-xIter)*(p_xTx0-1)+xdT*p_xdTx0;
        
        xIter = xIter - err/Jac;
        X0Iter = [xIter, 0, 0,  ydInit2(xIter, 0, aux)];
        iterNum = iterNum +1;
    end
    
    sol_S = ode113(@(t,x)odePcr3bp(t,x, aux),[0 DT], X0Iter, opts);
    sol_I.x = sol_S.x;
    sol_I.y =  zeros(size(sol_S.y));
    for i = 1:length(sol_S.x)
        t_epoch = sol_S.x(i);
        sol_I.y([1,2],i) = initial2synodic(sol_S.y([1,2],i),t_epoch);
%         sol_I.y([2,4],i) = initial2synodic(sol_S.y([1,3],i),t_epoch);
    end
    subplot(2,2,1)
    plot(sol_S.y(1,:),sol_S.y(2,:)); axis equal
    subplot(2,2,2)
    plot(sol_S.x,sol_S.y(1,:),sol_S.x,sol_S.y(2,:))
    grid on; legend('x','y')
    subplot(2,2,3)
    plot(sol_I.y(1,:),sol_I.y(2,:)); axis equal
    pause(0.0001)
    
    x0(ii,1) = xIter;
    dxdJ(ii-1,1) = (x0(ii,1)- x0(ii-1,1))/dJ;
    
    [J, iterNum, xIter, dxdJ(ii-1,1)]
    if abs(dxdJ(ii-1,1))<2e-1
        1;
    end
end

figure(101)
plot(Jmat(1:length(dxdJ)),dxdJ,'-o'); hold on; grid on
% end

