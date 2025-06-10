function [value, isterminal, direction] = event_impact(t, xx,aux,theta0)
%
% ײ������
if size(xx,2)~=1
    xx = xx';
end
if nargin == 3 || nargin == 4


    % ��ֹ������ײ����
    eq1 = norm(xx(1:3) - [-aux.mu;0;0]) - (aux.planet.req(3)+100)/aux.LU;
    
    eq2 = norm(xx(1:3) - [1-aux.mu;0;0]) - (aux.planet.req(10)+50)/aux.LU;
    
    value = eq1*eq2 ; % ײ����
    
    isterminal = 0; % �Ƿ�ִ�л�����ֹ����(����)
    
    direction = 0; % ������ֹ���������򸺻����ɸ�����
end
end

