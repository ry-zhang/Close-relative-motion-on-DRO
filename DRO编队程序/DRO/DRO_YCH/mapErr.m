function [y, tp, D] = mapErr(x, J,mu)
% mapping error function
% x - scalar


yd = ydFromJ_x(x, J,mu);
X0 = [x 0 0 yd];

DT = 2*pi*6;
opts = odeset('RelTol',1e-13,'AbsTol',1e-20, 'Events',@secYeq0Stop);

[T,X, Te, Xe, Ie] = ode113(@(t,x)pcr3bp(t,x,mu),[0 DT], X0, opts);

y =  Xe(end,1)-x;
tp = Te(end);

D = zeros(1, 1);

epsilon = 1e-8;
% numerical gradient
xPert = x +epsilon;
ydPert = ydFromJ_x(xPert, J,mu);
XPert = [xPert 0 0 ydPert];

[~,~, ~, XeTmp, ~] = ode113(@(t,x)pcr3bp(t,x,mu),[0 DT], XPert, opts);

dx = XeTmp(end, 1)-Xe(end,1);
D = dx/epsilon;

D = D - 1;





