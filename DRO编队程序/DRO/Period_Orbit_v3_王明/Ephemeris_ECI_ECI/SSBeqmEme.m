
function xxdot = SSBeqmEme(t , xx , aux)
%
% eme2000动力学方程(仅考虑月球和太阳引力场)
% 
% 输入：
% t:                                  [1x1]              儒略日(注意单位是s)
% xx:                               [6x1]             状态(km , km/s)
% aux:      
%       aux.earth.mu        [1x1]        地球引力场数(km^3/s^2)
%       aux.moon.mu       [1x1]        月球引力场数(km^3/s^2)
%       aux.sun.mu          [1x1]        太阳引力场数(km^3/s^2)
%
% 输出：
% xxdot:                          [6x1]              动力学方程(km , km/s^2)      
% 
% 作者: 张晨, 中科院空间应用工程与技术中心
% chenzhang@csu.ac.cn
% 2019/12/30
% ---------------------------------------------------------

% % 载入mu
% mu_earth = aux.earth.mu;
% mu_moon = aux.moon.mu;
% mu_sun = aux.sun.mu;

% 载入mu
mu_earth = aux.planet.mu(3);
mu_moon = aux.planet.mu(10);
mu_sun = aux.planet.mu(11);

% s -> day
jdate = t / 86400;

% 卫星位置矢量
r_sc = xx(1:3);

% 卫星位置矢量模
rmag_sc = sqrt( xx(1)^2 + xx(2)^2 + xx(3)^2 );

% 太阳位置矢量
% rv_sun = ephEme(jdate, 11, 3 , aux.DE430);
rv_sun = ephEme_mex(jdate, 11, 3 , aux.DE430);
r_sun = rv_sun(1:3);
rmag_sun = sqrt( r_sun(1)^2 + r_sun(2)^2 + r_sun(3)^2 );

% 月球位置矢量
% rv_moon = ephEme(jdate, 10, 3 , aux.DE430);
rv_moon = ephEme_mex(jdate, 10, 3 , aux.DE430);
r_moon = rv_moon(1:3);
rmag_moon = sqrt( r_moon(1)^2 + r_moon(2)^2 + r_moon(3)^2 );

% 太阳指向卫星位置矢量
r_sun2sc = r_sc - r_sun;
rmag_sun2sc = sqrt( r_sun2sc(1)^2 + r_sun2sc(2)^2 + r_sun2sc(3)^2 );

% 卫星指向月球位置矢量
r_moon2sc = r_sc - r_moon;
rmag_moon2sc = sqrt( r_moon2sc(1)^2 + r_moon2sc(2)^2 + r_moon2sc(3)^2 );

% --------------------- 计算加速度 ---------------------
% 加速度地球部分
acc_earth = - mu_earth * r_sc / rmag_sc^3;

% 加速度月球部分
% acc_moon = - mu_moon * (r_moon / rmag_moon^3 + r_moon2sc / rmag_moon2sc^3);
acc_moon = - mu_moon * (r_moon2sc / rmag_moon2sc^3);
% 加速度太阳部分
% acc_sun = - mu_sun * (r_sun / rmag_sun^3 + r_sun2sc / rmag_sun2sc^3);
acc_sun = - mu_sun * ( r_sun2sc / rmag_sun2sc^3);
% 行星加速度

posVel = jplEph_new(jdate,  3, aux.DE430);
acc_earth2 = posVel(:,3)/ 86400^2;                 % 换算成km/s
acc_planet = acc_earth + acc_moon + acc_sun-acc_earth2;

% 卫星加速度
acc_sc = [0 ; 0 ; 0];

% compute integration vector
xxdot = [ xx(4)
    xx(5)
    xx(6)
    acc_planet(1) + acc_sc(1);
    acc_planet(2) + acc_sc(2);
    acc_planet(3) + acc_sc(3)];

end
