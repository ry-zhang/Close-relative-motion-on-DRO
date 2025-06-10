
function crtbpMarkEM
% �����Ĺ���ϵ�µĵ����������
% 2015/4/10
% Copyright(C) Chen Zhang
% -------------------------------------------------------------

% leoFlag = 1; % ��leo
leoFlag = 0; % ��leo

% lloFlag = 1; % ��llo
lloFlag = 0; % ��llo

liFlag = 1; % ��ƽ����
% liFlag = 0; % ��ƽ����

% -------------------------------------------------------------
mu = 0.0121506683;
LU = 384405; % [km]

rSun = 6.955e5; % ̫���뾶
rEarth = 6378; % ����뾶
rMoon = 1737; % ����뾶

altiLeo = 200; % ����ͣ������߶�
altiLLO = 200; % ����Ŀ�����߶�

% ------------------------------------------------------------
% �����������
[xs , ys , zs] = sphere(30) ;
surf(-mu+xs*rEarth / LU,  ys*rEarth / LU, zs*rEarth / LU, 'FaceColor', 'w') ; hold on;
surf(1-mu+xs*rMoon / LU, ys*rMoon / LU, zs*rMoon / LU, 'FaceColor', 'w') ; hold on;

text(-mu,0,0, 'Earth');
text(1-mu,0,0, 'Moon');

% ��LEO
if leoFlag == 1
    ri = (rEarth + altiLeo) / LU; % ����ͣ��������ľ�
    angle = 0 : 0.05 : 2*pi + 0.05 ;
    plot(-mu + ri*cos(angle), ri*sin(angle), 'k-.') ; hold on;
end

% ��LLO
if lloFlag == 1
    rf = (rMoon + altiLLO) / LU; % ����Ŀ�������ľ�
    plot(1-mu + rf*cos(angle), rf*sin(angle), 'k-.') ; hold on;
end

% ��ƽ����
if liFlag == 1
    [Li_pos , ~] = crtbpLi(mu);
    plot(Li_pos(1,1) , Li_pos(1,2) , 'k.' , 'linewidth' , 1);
    plot(Li_pos(2,1) , Li_pos(2,2) , 'k.' , 'linewidth' , 1);
    plot(Li_pos(3,1) , Li_pos(3,2) , 'k.' , 'linewidth' , 1);
    plot(Li_pos(4,1) , Li_pos(4,2) , 'k.' , 'linewidth' , 1);
    plot(Li_pos(5,1) , Li_pos(5,2) , 'k.' , 'linewidth' , 1);
end

axis equal;
grid on;
view(0 , 90);

end
