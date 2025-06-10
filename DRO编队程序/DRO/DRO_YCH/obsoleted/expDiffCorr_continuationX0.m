clear; clc;
% find periodic orbit

global mu
mu = 0.01211;

L = lagrange_points(mu);
r1 = sqrt((L(:,1)+mu).^2+L(:,2).^2);
r2 = sqrt((L(:,1)-1+mu).^2+L(:,2).^2);
JLp=  (L(:,1).^2+L(:,2).^2)+2*((1-mu)./r1+mu./r2);

%%
J0 = 2.99;
x0 = 0.83; 

%% 
xiter = x0;
J = J0;

x0 = 0:0.1:5;
   
 
iter= 1;

[y, tp, D] = mapErr(xiter, J);
while norm(y)>1e-6 && iter <1000
    xiter = (xiter- y/D);
    [y, tp, D] = mapErr(xiter, J);
    iter = iter +1;
end
xPdc = xiter;

ydPdc = ydFromJ_x(xPdc, J);




%%
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
DT = 100;
pert = 0;
[T,X] = ode113(@(t,x)pcr3bp(t,x),[0 DT], [xPdc+pert, 0, 0, ydPdc], opts);

r1 = sqrt((X(:,1)+mu).^2+X(:,2).^2);
r2 = sqrt((X(:,1)-1+mu).^2+X(:,2).^2);   
Jode =  (X(:,1).^2+X(:,2).^2)+2*((1-mu)./r1+mu./r2)-(X(:,3).^2+X(:,4).^2);

%%


plot(X(:,1), X(:,2)); hold on
plot(L(:,1), L(:,2), 's')
grid on
axis equal
plot(-mu,0,'o');
plot(1-mu,0,'o');




