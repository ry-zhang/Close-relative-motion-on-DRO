function [y, tp, D] = mapErr_z(X, aux)
% mapping error function
% x - scalar
% global mu


X0 = X;

DT = 20;
opts = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events',@secYeq0Stop);

[T,~, Te, Xe, Ie] = ode113(@(t,x)odeR3bpAug(t,x,aux),[0 DT], X0, opts);

y =  Xe(end,:)-X0; %
tp = Te(end);

epsilon = 1e-8;
% numerical gradient
% 对x的数值求导
xPert = X(1) + epsilon;
ydPert = ydInit2(xPert, X(3), 0, aux);
XPert = [xPert 0 X(3) 0 ydPert 0];
[~,~, ~, XeTmp, ~] = ode113(@(t,x)odeR3bpAug(t,x,aux),[0 DT], XPert, opts);
dXdx = XeTmp(end, :)-XPert - y;

% 对z的数值求导
zPert = X(3) + epsilon;
ydPert = ydInit2(X(1), zPert, 0, aux);
XPert = [X(1) 0 zPert 0 ydPert 0];
[~,~, ~, XeTmp, ~] = ode113(@(t,x)odeR3bpAug(t,x,aux),[0 DT], XPert, opts);
dXdz = XeTmp(end, :)-XPert - y;

dXdxz = [dXdx;dXdz];
D = dXdxz/epsilon;

y = y([1:6]);
D = D(:,1:6);





