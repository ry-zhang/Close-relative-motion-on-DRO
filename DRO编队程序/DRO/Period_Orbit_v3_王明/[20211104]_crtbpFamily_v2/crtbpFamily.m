function crtbpFamily
%
% 绘制三体周期轨道族
%
% 作者：张晨，王明，张皓
% 单位：中科院空间应用工程与技术中心
% 时间：2021年3月30日
% chenzhang@csu.ac.cn
% -----------------------------------------------------------

% LOP-G飞行器

clc; clear; close all;

% 地月引力常数
mu = 0.01215058560962404;

% 地月LU
LU = 384405;

% 地月TU
TU = 375676.96752;

% ------------------------- DRO --------------------------
% FileName = 'DRO_2D.txt';
% FileName = 'DRO_3D_North.txt';
% FileName = 'DRO_3D_South.txt';

% ------------------------- DPO --------------------------
% FileName = 'DPO_2D.txt';
% FileName = 'DPO_3D_North.txt';
% FileName = 'DPO_3D_South.txt';

% ------------------------ LoPO --------------------------
% FileName = 'LoPO_2D.txt';

% ---------------------------- L1 --------------------------
% FileName = 'L1_Lya_2D.txt';
% FileName = 'L1_Halo_3D_North.txt';
% FileName = 'L1_Halo_3D_South.txt';

% ---------------------------- L2 --------------------------
% FileName = 'L2_Lya_2D.txt';
% FileName = 'L2_Halo_3D_North.txt';
% FileName = 'L2_Halo_3D_South.txt';

% --------------------------- L3 ---------------------------
% FileName = 'L3_Lya_2D.txt';
% FileName = 'L3_Halo_3D_North.txt';
% FileName = 'L3_Halo_3D_South.txt';

% --------------------------- L4 ---------------------------

% FileName = 'L4_short_2D.txt';
% FileName = 'L4_Long_2D.txt';
% FileName = 'L4_Long_3D_North.txt';
% FileName = 'L4_Long_3D_South.txt';

% --------------------------- L5 ---------------------------
% FileName = 'L5_short_2D.txt';
% FileName = 'L5_Long_2D.txt';
% FileName = 'L5_Long_3D_North.txt';
% FileName = 'L5_Long_3D_South.txt';

% --------------------- Resonant -------------------------
% FileName = 'RO_1_1_2D.txt';
% FileName = 'RO_1_2_2D.txt';
% FileName = 'RO_1_3_2D.txt';
% FileName = 'RO_1_4_2D.txt';
% FileName = 'RO_2_1_2D.txt';
% FileName = 'RO_2_3_2D.txt';
% FileName = 'RO_3_1_2D.txt';
% FileName = 'RO_3_2_2D.txt';
FileName = 'RO_4_1_2D.txt';
% FileName = 'RO_4_3_2D.txt';

% ----------------------- 读取数据 -----------------------
fid = fopen(FileName , 'r');
tline = fgets(fid);
tline = fgets(fid);
count = 1;
while tline ~= -1
    xxNow = textscan(tline,'%f');
    xxNow = cell2mat(xxNow);
    X0Mtx(count,:) = xxNow';
    tline = fgets(fid);
    count = count +1 ;
end
fclose('all');
% X0Mtx = X0Mtx(100:200,:);
% min(X0Mtx(: , 7))
% max(X0Mtx(: , 7))

% min(X0Mtx(: , 8))
% max(X0Mtx(: , 8))

% 旋转系画图
h1 = figure(1); hold on; grid on; axis equal;
set(h1 , 'position' , [100 , 100 , 600 , 600]);

% 画EM系统
crtbpMarkEM;

% 该族轨道周期范围
Pmin = min(X0Mtx(: , end - 1));
Pmax = max(X0Mtx(: , end - 1));
fprintf('轨道周期：[%0.2f , %0.2f] days \n' , Pmin  , Pmax )
% fprintf('轨道周期：[%0.2f , %0.2f] days \n' , Pmin * TU / 86400 , Pmax * TU / 86400)

% 该族轨道能量范围
Cmin = min(X0Mtx(: , end));
Cmax = max(X0Mtx(: , end));

fprintf('轨道能量：[%0.2f , %0.2f] \n' , Cmin , Cmax)

for iLoop = 1 : size(X0Mtx , 1)
    
    % 周期轨道初值
    xx_periOrb = X0Mtx(iLoop , 1:6);
    
    % 周期轨道周期
    P_periOrb = X0Mtx(iLoop , 7);
    
    % 周期轨道能量
    C_periOrb = X0Mtx(iLoop , 8);
    
    % 画周期轨道
    crtbpFlow3D(xx_periOrb , P_periOrb , C_periOrb , Cmin , Cmax , mu);
    
end

view(45 , 45);

end

function [tt , xx] = crtbpFlow3D(xx0 , P , C , Cmin , Cmax , mu)
%
% crtbp三维数值积分
%
% 例子：
% crtbpMarkEM
% aux = crtbpAuxEM;
% xx0 = [0.75 , 0 , 0 , 0.15 , 0 , 0];
% crtbpFlow3D(xx0 , 10 , aux.mu)
%
% 轨道颜色
% jet / hsv / parula
% 
% 作者：张晨
% 单位：中科院空间应用工程与技术中心
% 时间：2021年3月30日
% chenzhang@csu.ac.cn
% -----------------------------------------------------------

% 颜色插值，从[Cmin , Cmax]区间映射到[1 , 256]区间，得到对应RGB颜色
Nmin = 1;
Nmax = 256;
colorTemp = interp1([1 : 256]' , parula , (Nmax - Nmin) / (Cmax - Cmin) * (C - Cmin) + Nmin , 'spline');

% 颜色超过[0,1]边界，修正
colorTemp(colorTemp > 1) = 1;
colorTemp(colorTemp < 0) = 0;

% 数值积分，画轨道
options = odeset('Reltol', 1e-12, 'AbsTol',1e-12);
[tt , xx] = ode113(@fx , [0 , P] , xx0 , options , mu);
plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'color' , colorTemp , 'linewidth' , 1);

end


function dxxdt = fx (t, xx, mu)
%
% 圆型限制性三体问题（CRTBP）动力学方程
%
% 作者：张晨
% 单位：中科院空间应用工程与技术中心
% 时间：2021年3月30日
% chenzhang@csu.ac.cn
% -----------------------------------------------------------
r1cube = ((xx(1) + mu)^2 + xx(2)^2 + xx(3)^2)^(1.5);
r2cube = ((xx(1) - 1 + mu)^2 + xx(2)^2 + xx(3)^2)^(1.5);
dxxdt = [xx(4);
    xx(5);
    xx(6);
    xx(1) + 2*xx(5) - (1 - mu) * (xx(1) + mu) / r1cube - mu*(xx(1) - 1 + mu) / r2cube;
    xx(2) - 2*xx(4) - (1 - mu) * xx(2) / r1cube - mu * xx(2) / r2cube;
    -(1 - mu) * xx(3) / r1cube - mu * xx(3) / r2cube];
end

function crtbpMarkEM
%
% 画地球和月球
%
% 作者：张晨
% 单位：中科院空间应用工程与技术中心
% 时间：2021年3月30日
% chenzhang@csu.ac.cn
% -----------------------------------------------------------

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

text( -mu,  0, 0, 'Earth') ; hold on;
text(1 - mu,  0, 0, 'Moon') ; hold on;

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
    plot(Li_pos(1,1) , Li_pos(1,2) , 'r+' , 'linewidth' , 1 , 'markersize' , 5);
    plot(Li_pos(2,1) , Li_pos(2,2) , 'r+' , 'linewidth' , 1 , 'markersize' , 5);
    plot(Li_pos(3,1) , Li_pos(3,2) , 'r+' , 'linewidth' , 1 , 'markersize' , 5);
    plot(Li_pos(4,1) , Li_pos(4,2) , 'r+' , 'linewidth' , 1 , 'markersize' , 5);
    plot(Li_pos(5,1) , Li_pos(5,2) , 'r+' , 'linewidth' , 1 , 'markersize' , 5);
    
    text(Li_pos(1,1) , Li_pos(1,2) , 'L1' , 'linewidth' , 1);
    text(Li_pos(2,1) , Li_pos(2,2) , 'L2' , 'linewidth' , 1);
    text(Li_pos(3,1) , Li_pos(3,2) , 'L3' , 'linewidth' , 1);
    text(Li_pos(4,1) , Li_pos(4,2) , 'L4' , 'linewidth' , 1);
    text(Li_pos(5,1) , Li_pos(5,2) , 'L5' , 'linewidth' , 1);
end

axis equal;
grid on;
xlabel('x (LU)');
ylabel('y (LU)');
zlabel('z (LU)');
view(0 , 90);

end

function [Li_pos , Li_c] = crtbpLi(mu)
% 计算crtbp的5个拉格朗日点和能量
%
% 例子：
% % [1]
% aux = crtbpAuxEM;
% [Li_pos , Li_c] = crtbpLi(aux.mu)
%
% % [2]
% aux = crtbpAuxSE;
% [Li_pos , Li_c] = crtbpLi(aux.mu)
%
% 作者：张晨
% 单位：中科院空间应用工程与技术中心
% 时间：2021年3月30日
% chenzhang@csu.ac.cn
% -----------------------------------------------------------
% 设置初始猜测值
x012 = (mu/3)^(1/3) ;
x03  = 1-(7/12)*mu  ;
% 计算5个平动点坐标
options = optimset('TolFun', 2.5e-14, 'TolX', 2.5e-14) ;
L1  = 1 - mu - fzero(@(x) x^5-(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2+2*mu*x-mu, ...
    x012, options) ;
L2  = 1 - mu + fzero(@(x) x^5+(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2-2*mu*x-mu, ...
    x012, options) ;
L3  =   - mu - fzero(@(x) x^5+(2+mu)*x^4+(1+2*mu)*x^3-(1-mu)*x^2-2*(1-mu)*x...
    -(1-mu), x03,  options) ;
L4x = 0.5 - mu ;
L4y = sqrt(3)/2 ;
% 构造5个平动点位置
Li_pos = [L1, 0;
    L2, 0;
    L3, 0;
    L4x, L4y;
    L4x, -L4y];
% 如果输出大于1则计算5个平动点能量
if nargout > 1
    xx = [Li_pos , zeros(5 , 2)];
    Li_c = crtbpJacobi2D(xx , mu);
end
end

function CC = crtbpJacobi2D (xx, mu)
% 计算雅可比能量C
%
% 输入：
% xx       [nx4]
% 输出：
% C        [nx1]
%
% 作者：张晨
% 单位：中科院空间应用工程与技术中心
% 时间：2021年3月30日
% chenzhang@csu.ac.cn
% -----------------------------------------------------------
if size(xx , 2) ~= 4
    disp('xx是n * 4')
end
r1 = sqrt((xx(: , 1)+mu).^2   + xx(: , 2).^2) ;
r2 = sqrt((xx(: , 1)+mu-1).^2 + xx(: , 2).^2) ;
CC = -(xx(: , 3).^2 + xx(: , 4).^2) ...
    + (xx(: , 1).^2 + xx(: , 2).^2) ...
    + 2 * (1 - mu) ./ r1 + 2 * mu ./ r2;
end
