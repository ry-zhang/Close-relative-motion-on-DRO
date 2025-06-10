function [y, tp, D] = mapErr_xz(x, z, aux)
% mapping error function
% x - scalar
% global mu

xd = 0;
yd = ydInit2(x, z, xd, aux);

X0 = [x 0 z 0 yd 0];

DT = 2*pi*3;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events',@secYeq0Stop);

[T,X, Te, Xe, Ie] = ode113(@(t,x)odeR3bpAug(t,x,aux),[0 DT], X0, opts);

y =  Xe(end,1)-x; %
tp = Te(end);

D = zeros(1, 1);

epsilon = 1e-8;
% numerical gradient
xPert = x + epsilon;
ydPert = ydInit2(xPert, z, 0, aux);
XPert = [xPert 0 z 0 ydPert 0];

[~,~, ~, XeTmp, ~] = ode113(@(t,x)odeR3bpAug(t,x, aux),[0 DT], XPert, opts);

dy = XeTmp(end, 1)-Xe(end,1);
D = dy/epsilon;

D = D - 1;





