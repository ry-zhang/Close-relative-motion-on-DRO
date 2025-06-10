function [value, isterminal, direction] = event_impact(t, xx,aux,theta0)
%
% 撞地球检测
if size(xx,2)~=1
    xx = xx';
end
if nargin == 3 || nargin == 4


    % 终止条件（撞地球）
    eq1 = norm(xx(1:3) - [-aux.mu;0;0]) - (aux.planet.req(3)+100)/aux.LU;
    
    eq2 = norm(xx(1:3) - [1-aux.mu;0;0]) - (aux.planet.req(10)+50)/aux.LU;
    
    value = eq1*eq2 ; % 撞地球
    
    isterminal = 0; % 是否执行积分终止条件(单次)
    
    direction = 0; % 积分终止方向，由正向负或者由负向正
end
end

