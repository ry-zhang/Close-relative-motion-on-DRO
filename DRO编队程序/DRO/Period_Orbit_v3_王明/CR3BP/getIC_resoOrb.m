
function [yy , aux] = getIC_resoOrb(aux)
%
% ����ಽ��й�������ֵ
%
% �ο���
% [1] Mar Vaquero, Leveraging Resonant-Orbit Manifolds to Design Transfers
% Between Libration-Point Orbits, 2014
%
% ����: �ų�, �п�Ժ�ռ�Ӧ�ù����뼼������
% chenzhang@csu.ac.cn
% 2021/06/07
% -----------------------------------------------------------

mu = aux.mu;
LU = aux.LU;
VU = aux.VU;
TU = aux.TU;

% ----------------------- ���ò��� -----------------------
% ����ڵ���
node_n = aux.node_n;

% ����Ȧ��
resoOrb_p = aux.resoOrb_p/aux.gcd;

% С����Ȧ��
resoOrb_q = aux.resoOrb_q/aux.gcd;

% С��������
Period_q = 2 * pi;

% ���ص�߶�
resoOrb_rp = aux.resoOrb_rp;

% ���ǹ�����
resoOrb_incl = aux.resoOrb_incl;

% ����������ྭ
resoOrb_raan = aux.resoOrb_raan;

% ����������
resoOrb_tanom = aux.resoOrb_tanom;

% ------------------------- �����ֵ ---------------------------
% ���������� (TU)
Period_p = Period_q * resoOrb_q / resoOrb_p;

% �볤�� (LU)
sma_res = ((1 - mu) * (Period_p / (2 * pi))^2) ^ (1/3);

% ���ع�������
beta_res = sqrt(2 * (1 - mu) / resoOrb_rp - (1 - mu) / sma_res) / ...
    sqrt((1 - mu) / resoOrb_rp);

% �ܻ���ʱ��
P_periOrb = resoOrb_q * Period_q;

% ����ϵ��ֵ
sI = sin(resoOrb_incl);
cI = cos(resoOrb_incl);
sR = sin(resoOrb_raan);
cR = cos(resoOrb_raan);
sT = sin(resoOrb_tanom);
cT = cos(resoOrb_tanom);
xi_ini = [resoOrb_rp * (cT * cR - sT * sR * cI);
    resoOrb_rp * (cT * sR + sT * cR * cI);
    resoOrb_rp * (sT * sI);
    beta_res * -sqrt((1 - mu) / resoOrb_rp) * (sT * cR + cT * sR * cI);
    beta_res * -sqrt((1 - mu) / resoOrb_rp) * (sT * sR - cT * cR * cI) ;
    beta_res * sqrt((1 - mu) / resoOrb_rp) * cT * sI];

% ------------------------- ���� & ��ֵ ---------------------------
% ����ϵ����
options = odeset('Reltol' , aux.tol , 'AbsTol' , aux.tol);
[tt_ini , xx_ini] = ode113(@eqm2b , [0 , P_periOrb] , xi_ini , options , 1 - mu);
xx_ini_interp = interp1(tt_ini , xx_ini , linspace(0 , P_periOrb , node_n + 1) , 'spline');
% ����ϵ -> ��תϵ
[xx_rot , tt_rot] = crtbpP1centered2synodic3D(xx_ini , tt_ini , mu);
if aux.UseSymmetric_IO == 1
    % ��ֵ����
    aux.periOrb_P = P_periOrb/2;
    options = odeset('Reltol', 1e-12 , 'AbsTol' , 1e-12);
    x0 = xx_rot(1,:);
    [tt_rot , xx_rot] = ode113(@crtbpEqm3D, [0 , aux.periOrb_P] , x0 , options , aux);
    xx_rot_interp = interp1(tt_rot , xx_rot , linspace(0 , aux.periOrb_P , node_n + 1) , 'spline');
    temp = xx_rot_interp(1 : node_n , :)';
    yy = [temp(:) ; aux.periOrb_P];
    yy([2,4,6]) = [];                                               % ������öԳ��ԣ����ɱ���ȥ��y, x_dot, z_dot
else
    % ��ֵ
    aux.periOrb_P = P_periOrb;
    xx_rot_interp = interp1(tt_rot , xx_rot , linspace(0 , aux.periOrb_P , node_n + 1) , 'spline');
    
    % ����ಽ��г�ֵ
    temp = xx_rot_interp(1 : node_n , :)';
    yy = [temp(:) ; aux.periOrb_P];
    
    % �̶���������һ�����yֵ��
    aux.yFix = yy(2);
end



% ----------------------- ��ͼ -----------------------
h1 = figure(1); hold on; grid on; axis equal;
set(h1 , 'position' , [100 , 100 , 600 , 400]);
plot_o([0 , 0 , 0] , 1 , 'k--');
text(0,0,0 , 'Earth');
text(1,0,0 , 'Moon');

[xs , ys , zs] = sphere(30) ;
surf(6378 * xs / LU,  6378 * ys / LU, 6378 * zs / LU, 'FaceColor', 'w') ; hold on;
surf(1738 * xs / LU + 1 , 1738 * ys / LU, 1738 * zs / LU, 'FaceColor', 'w') ; hold on;

aux.h1_g = plot3(xx_ini(: , 1) , xx_ini(: , 2) , xx_ini(: , 3) , 'g' , 'linewidth' , 1);
plot3(xx_ini(1 , 1) , xx_ini(1 , 2) , xx_ini(1 , 3) , 'g.' , 'linewidth' , 3 , 'markersize' , 20);
plot3(xx_ini_interp(: , 1) , xx_ini_interp(: , 2) , xx_ini_interp(: , 3) , 'go' , 'linewidth' , 1);

xlabel('x/LU');
ylabel('y/LU');
zlabel('z/LU');
title('ini frame')

% ---------------------- ��תϵ��ͼ ---------------------
h2 = figure(2); hold on; grid on;
set(h2 , 'position' , [200 , 200 , 600 , 400]);

crtbpMarkEM;
axis equal;

aux.h2_g = plot3(xx_rot(: , 1) , xx_rot(: , 2) , xx_rot(: , 3) , 'g' , 'linewidth' , 1);
plot3(xx_rot(1 , 1) , xx_rot(1 , 2) , xx_rot(1 , 3) , 'g.' , 'linewidth' , 3 , 'markersize' , 20);
plot3(xx_rot_interp(: , 1) , xx_rot_interp(: , 2) , xx_rot_interp(: , 3) , 'go' , 'linewidth' , 1);

xlabel('x/LU');
ylabel('y/LU');
zlabel('z/LU');
title('rot frame')


end


function xxdot = eqm2b (t , xx , mu)
% ���Ƕ��嶯��ѧ����
%
% ���룺
% t:              [1x1]             ʱ��(s)
% xx:             [6x1]             ״̬(km , km/s)
% mu:         [1x1]              ������������ (km^3/s^2)
%
% �����
% xxdot:         [6x1]            ����ѧ����(km , km/s^2)
%
% ����: �ų�, �п�Ժ�ռ�Ӧ�ù����뼼������
% chenzhang@csu.ac.cn
% 2019/12/30
% -----------------------------------------------------------

r2 = xx(1) * xx(1) + xx(2) * xx(2) + xx(3) * xx(3);
r1 = sqrt(r2);
r3 = r2 * r1;
xxdot = [xx(4);
    xx(5);
    xx(6);
    xx(1) * ( - mu / r3);
    xx(2) * ( - mu / r3);
    xx(3) * ( - mu / r3)];

end
