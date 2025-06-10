
function xxdot = SSBeqmEme(t , xx , aux)
%
% eme2000����ѧ����(�����������̫��������)
% 
% ���룺
% t:                                  [1x1]              ������(ע�ⵥλ��s)
% xx:                               [6x1]             ״̬(km , km/s)
% aux:      
%       aux.earth.mu        [1x1]        ������������(km^3/s^2)
%       aux.moon.mu       [1x1]        ������������(km^3/s^2)
%       aux.sun.mu          [1x1]        ̫����������(km^3/s^2)
%
% �����
% xxdot:                          [6x1]              ����ѧ����(km , km/s^2)      
% 
% ����: �ų�, �п�Ժ�ռ�Ӧ�ù����뼼������
% chenzhang@csu.ac.cn
% 2019/12/30
% ---------------------------------------------------------

% % ����mu
% mu_earth = aux.earth.mu;
% mu_moon = aux.moon.mu;
% mu_sun = aux.sun.mu;

% ����mu
mu_earth = aux.planet.mu(3);
mu_moon = aux.planet.mu(10);
mu_sun = aux.planet.mu(11);

% s -> day
jdate = t / 86400;

% ����λ��ʸ��
r_sc = xx(1:3);

% ����λ��ʸ��ģ
rmag_sc = sqrt( xx(1)^2 + xx(2)^2 + xx(3)^2 );

% ̫��λ��ʸ��
% rv_sun = ephEme(jdate, 11, 3 , aux.DE430);
rv_sun = ephEme_mex(jdate, 11, 3 , aux.DE430);
r_sun = rv_sun(1:3);
rmag_sun = sqrt( r_sun(1)^2 + r_sun(2)^2 + r_sun(3)^2 );

% ����λ��ʸ��
% rv_moon = ephEme(jdate, 10, 3 , aux.DE430);
rv_moon = ephEme_mex(jdate, 10, 3 , aux.DE430);
r_moon = rv_moon(1:3);
rmag_moon = sqrt( r_moon(1)^2 + r_moon(2)^2 + r_moon(3)^2 );

% ̫��ָ������λ��ʸ��
r_sun2sc = r_sc - r_sun;
rmag_sun2sc = sqrt( r_sun2sc(1)^2 + r_sun2sc(2)^2 + r_sun2sc(3)^2 );

% ����ָ������λ��ʸ��
r_moon2sc = r_sc - r_moon;
rmag_moon2sc = sqrt( r_moon2sc(1)^2 + r_moon2sc(2)^2 + r_moon2sc(3)^2 );

% --------------------- ������ٶ� ---------------------
% ���ٶȵ��򲿷�
acc_earth = - mu_earth * r_sc / rmag_sc^3;

% ���ٶ����򲿷�
% acc_moon = - mu_moon * (r_moon / rmag_moon^3 + r_moon2sc / rmag_moon2sc^3);
acc_moon = - mu_moon * (r_moon2sc / rmag_moon2sc^3);
% ���ٶ�̫������
% acc_sun = - mu_sun * (r_sun / rmag_sun^3 + r_sun2sc / rmag_sun2sc^3);
acc_sun = - mu_sun * ( r_sun2sc / rmag_sun2sc^3);
% ���Ǽ��ٶ�

posVel = jplEph_new(jdate,  3, aux.DE430);
acc_earth2 = posVel(:,3)/ 86400^2;                 % �����km/s
acc_planet = acc_earth + acc_moon + acc_sun-acc_earth2;

% ���Ǽ��ٶ�
acc_sc = [0 ; 0 ; 0];

% compute integration vector
xxdot = [ xx(4)
    xx(5)
    xx(6)
    acc_planet(1) + acc_sc(1);
    acc_planet(2) + acc_sc(2);
    acc_planet(3) + acc_sc(3)];

end
