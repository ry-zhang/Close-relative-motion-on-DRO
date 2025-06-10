
function dxxdt = crtbpEqm3D(t , xx , aux)
% crtbp����ѧ����
%
% ���룺
% t:               [1x1]            ʱ��(s)
% xx:             [6x1]            ״̬(km , km/s)
% mu:           [1x1]            �������� (km^3/s^2)
%% �����
% xxdot:        [6x1]            ����ѧ����(km , km/s^2)
%
% ����
% xx0 = [0.75 , 0 , 0 , 0.15 , 0 , 0];
% options = odeset('RelTol',1e-8,'AbsTol', 1e-8);
% [tt , xx]=ode113(@crtbpEqm3D , [0 , 10] , xx0 , options , aux);
%
% ����: �ų�, �п�Ժ�ռ�Ӧ�ù����뼼������
% chenzhang@csu.ac.cn
% 2019/12/30
% ---------------------------------------------------------
mu = aux.mu;

r1cube = ((xx(1) + mu)^2 + xx(2)^2 + xx(3)^2)^(1.5);
r2cube = ((xx(1) - 1 + mu)^2 + xx(2)^2 + xx(3)^2)^(1.5);
dxxdt = [xx(4);
    xx(5);
    xx(6);
    xx(1) + 2*xx(5) - (1 - mu) * (xx(1) + mu) / r1cube - mu*(xx(1) - 1 + mu) / r2cube;
    xx(2) - 2*xx(4) - (1 - mu) * xx(2) / r1cube - mu * xx(2) / r2cube;
    -(1 - mu) * xx(3) / r1cube - mu * xx(3) / r2cube];

end

