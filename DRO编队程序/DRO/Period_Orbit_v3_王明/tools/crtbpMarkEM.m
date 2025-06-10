
function crtbpMarkEM
% 画地心惯性系下的地球和月球轨道
% 2015/4/10
% Copyright(C) Chen Zhang
% -------------------------------------------------------------

% leoFlag = 1; % 画leo
leoFlag = 0; % 画leo

% lloFlag = 1; % 画llo
lloFlag = 0; % 画llo

liFlag = 1; % 画平动点
% liFlag = 0; % 画平动点

% -------------------------------------------------------------
mu = 0.0121506683;
LU = 384405; % [km]

rSun = 6.955e5; % 太阳半径
rEarth = 6378; % 地球半径
rMoon = 1737; % 月球半径

altiLeo = 200; % 近地停泊轨道高度
altiLLO = 200; % 近月目标轨道高度

% ------------------------------------------------------------
% 画地球和月球
[xs , ys , zs] = sphere(30) ;
surf(-mu+xs*rEarth / LU,  ys*rEarth / LU, zs*rEarth / LU, 'FaceColor', 'w') ; hold on;
surf(1-mu+xs*rMoon / LU, ys*rMoon / LU, zs*rMoon / LU, 'FaceColor', 'w') ; hold on;

text(-mu,0,0, 'Earth');
text(1-mu,0,0, 'Moon');

% 画LEO
if leoFlag == 1
    ri = (rEarth + altiLeo) / LU; % 近地停泊轨道地心距
    angle = 0 : 0.05 : 2*pi + 0.05 ;
    plot(-mu + ri*cos(angle), ri*sin(angle), 'k-.') ; hold on;
end

% 画LLO
if lloFlag == 1
    rf = (rMoon + altiLLO) / LU; % 近月目标轨道月心距
    plot(1-mu + rf*cos(angle), rf*sin(angle), 'k-.') ; hold on;
end

% 画平动点
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
