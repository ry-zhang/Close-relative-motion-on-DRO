function aux = initialize(aux)
%-------------初始化星历参数---------------------------------------------------
%赤道系/黄道系转换矩阵
rotation_matrices;

% 天体引力参数from DE430
xmu(1) = 2.2031780000000021E+04 ; % mercury
xmu(2) = 3.2485859200000006E+05;
xmu(3) = 398600.435436;
xmu(4) = 4.2828375214000022E+04; %mars
xmu(5) = 1.2671276480000021E+08 ; % jupiter
xmu(6) =  3.7940585200000003E+07;   % saturn
xmu(7) = 5.7945486000000080E+06 ;
xmu(8) = 6.8365271005800236E+06 ;
xmu(9) = 9.7700000000000068E+02; %pluto
xmu(10) = 4902.801;  %moon
xmu(11) = 1.3271244004193938E+11; % sun, km^3*s^-2

% 日地月系统参数
aux.LU = 384000;% km
aux.TU = 3.746047897699050e+05;
aux.VU = 1.025080325950626;
aux.mu = 0.01215;

% 需要配置天体引力参数/开始儒略日/坐标转换矩阵/考虑的天体类型
aux.xmu = xmu;
aux.ICRF_2_MeanEclpJ2k = ICRF_2_MeanEclpJ2k;
aux.ICRF_2_J2K = ICRF_2_J2K;
% aux.threeBodyList = [1:2, 4:6]; %考虑的三体类型： 水 金 火 土 木星
aux.threeBodyList = [];% 只考虑地月日

% % 初始历元
% aux.t0UTC  = [2030 1 1 12 0 0];%初始时刻
% 从UTC计算JD
UTC2TCB  = [0,0,0,0, 1, 9.186]; % @2020+之后的跳秒

% 转儒略日
t0TCB = aux.t0UTC+UTC2TCB;
jd0 = juliandate(t0TCB(1), t0TCB(2), t0TCB(3), t0TCB(4), t0TCB(5), t0TCB(6));%ICRF坐标系,day
aux.jd0 = jd0;% day
end