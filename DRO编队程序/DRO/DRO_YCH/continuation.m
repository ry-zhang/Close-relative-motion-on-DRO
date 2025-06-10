%% cost function
function [delta_T,J_period,state_ini] = continuation(J,x0,T_DRO_desired,opts1,mu)

iter = 1;
xiter = x0;
[y, ~, D] = mapErr(xiter, J,mu);
while norm(y)>1e-13 && iter <50
    xiter = (xiter - y/D);
    [y, ~, D] = mapErr(xiter, J,mu);
    iter = iter +1;
end

xPdc = xiter;
ydPdc = ydFromJ_x(xPdc, J,mu);
DT = 6*2*pi;
sol = ode113(@(t,x)pcr3bp(t,x,mu),[0 DT], [xPdc, 0, 0, ydPdc], opts1);
T_temp = diff(sol.xe);

if ~isempty(T_temp)
    T0 = T_temp(1);
else
    T0 = T_DRO_desired/2;
end

delta_T = T0/(2*pi)-T_DRO_desired; % cost function
delta_T = 1e3*delta_T;

if max(T_temp-mean(T_temp))>1e-4
    delta_T = delta_T + max(T_temp-mean(T_temp));
end

J_period = [J,T0/(2*pi)];
state_ini = [xPdc, 0, 0, ydPdc];
end