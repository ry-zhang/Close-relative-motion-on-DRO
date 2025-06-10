
function [yy , aux] = getIC_resoOrb(aux)
%
% 计算多步打靶共振轨道初值
%
% 参考：
% [1] Mar Vaquero, Leveraging Resonant-Orbit Manifolds to Design Transfers
% Between Libration-Point Orbits, 2014
%
% 作者: 张晨, 中科院空间应用工程与技术中心
% chenzhang@csu.ac.cn
% 2021/06/07
% -----------------------------------------------------------

mu = aux.mu;
LU = aux.LU;
VU = aux.VU;
TU = aux.TU;

% ----------------------- 配置参数 -----------------------
% 轨道节点数
node_n = aux.node_n;

% 卫星圈数
resoOrb_p = aux.resoOrb_p/aux.gcd;

% 小天体圈数
resoOrb_q = aux.resoOrb_q/aux.gcd;

% 小天体周期
Period_q = 2 * pi;

% 近地点高度
resoOrb_rp = aux.resoOrb_rp;

% 卫星轨道倾角
resoOrb_incl = aux.resoOrb_incl;

% 卫星升交点赤经
resoOrb_raan = aux.resoOrb_raan;

% 卫星真近点角
resoOrb_tanom = aux.resoOrb_tanom;

% ------------------------- 计算初值 ---------------------------
% 共振轨道周期 (TU)
Period_p = Period_q * resoOrb_q / resoOrb_p;

% 半长轴 (LU)
sma_res = ((1 - mu) * (Period_p / (2 * pi))^2) ^ (1/3);

% 近地轨道脉冲比
beta_res = sqrt(2 * (1 - mu) / resoOrb_rp - (1 - mu) / sma_res) / ...
    sqrt((1 - mu) / resoOrb_rp);

% 总积分时间
P_periOrb = resoOrb_q * Period_q;

% 惯性系初值
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

% ------------------------- 积分 & 插值 ---------------------------
% 惯性系积分
options = odeset('Reltol' , aux.tol , 'AbsTol' , aux.tol);
[tt_ini , xx_ini] = ode113(@eqm2b , [0 , P_periOrb] , xi_ini , options , 1 - mu);
xx_ini_interp = interp1(tt_ini , xx_ini , linspace(0 , P_periOrb , node_n + 1) , 'spline');
% 惯性系 -> 旋转系
[xx_rot , tt_rot] = crtbpP1centered2synodic3D(xx_ini , tt_ini , mu);
if aux.UseSymmetric_IO == 1
    % 数值积分
    aux.periOrb_P = P_periOrb/2;
    options = odeset('Reltol', 1e-12 , 'AbsTol' , 1e-12);
    x0 = xx_rot(1,:);
    [tt_rot , xx_rot] = ode113(@crtbpEqm3D, [0 , aux.periOrb_P] , x0 , options , aux);
    xx_rot_interp = interp1(tt_rot , xx_rot , linspace(0 , aux.periOrb_P , node_n + 1) , 'spline');
    temp = xx_rot_interp(1 : node_n , :)';
    yy = [temp(:) ; aux.periOrb_P];
    yy([2,4,6]) = [];                                               % 如果利用对称性，自由变量去掉y, x_dot, z_dot
else
    % 插值
    aux.periOrb_P = P_periOrb;
    xx_rot_interp = interp1(tt_rot , xx_rot , linspace(0 , aux.periOrb_P , node_n + 1) , 'spline');
    
    % 输出多步打靶初值
    temp = xx_rot_interp(1 : node_n , :)';
    yy = [temp(:) ; aux.periOrb_P];
    
    % 固定共振轨道第一个点的y值！
    aux.yFix = yy(2);
end



% ----------------------- 画图 -----------------------
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

% ---------------------- 旋转系画图 ---------------------
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
% 卫星二体动力学方程
%
% 输入：
% t:              [1x1]             时间(s)
% xx:             [6x1]             状态(km , km/s)
% mu:         [1x1]              地球引力常数 (km^3/s^2)
%
% 输出：
% xxdot:         [6x1]            动力学方程(km , km/s^2)
%
% 作者: 张晨, 中科院空间应用工程与技术中心
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
