function aux = initialize(aux)
%-------------初始化星历参数---------------------------------------------------

%赤道系/黄道系转换矩阵
rotation_matrices;

% 天体引力参数from DE430
xmu(1) = 2.2031780000000021E+04 ; % mercury
xmu(2) = 3.2485859200000006E+05; % venus
xmu(3) = 398600.435436; % earth-moon
xmu(4) = 4.2828375214000022E+04; % mars
xmu(5) = 1.2671276480000021E+08 ; % jupiter
xmu(6) =  3.7940585200000003E+07; % saturn
xmu(7) = 5.7945486000000080E+06 ; % Uranus
xmu(8) = 6.8365271005800236E+06 ; % Neptune
xmu(9) = 9.7700000000000068E+02; % pluto
xmu(10) = 4902.801;  % moon
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
% UTC2TDB  = [0,0,0,0, 1, 9.186]; % @2020+之后的跳秒,juliandate不需要跳秒
% t0TCB = aux.t0UTC+UTC2TCB;

% 转儒略日
isjd0 = isfield(aux,'jd0') && ~isempty(aux.jd0);
ist0UTC = isfield(aux,'t0UTC') && ~isempty(aux.t0UTC);
if isjd0 && ~ist0UTC
    aux.t0TDB = datetime(aux.jd0,'ConvertFrom','juliandate','Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
    aux.t0UTC = datetime(aux.t0TDB,'Format','yyyy-MM-dd HH:mm:ss.SSSSSS') - minutes(1) - seconds(9.186);
elseif ~isjd0 && ist0UTC
    aux.t0TDB = aux.t0UTC + minutes(1) + seconds(9.186);
    aux.jd0 = juliandate(aux.t0TDB); %ICRF坐标系,day
elseif isjd0 && ist0UTC
    aux.t0TDB = datetime(aux.jd0,'ConvertFrom','juliandate');
    t0TDB2 = aux.t0UTC + minutes(1) + seconds(9.186);
    if abs(aux.jd0 - juliandate(t0TDB2))>1e-8
        warning('jd0与t0UTC不相符，请注意jd0为质心动力学时(Barycentric Dynamical Time),在此以jd0为准')
    end
elseif ~isjd0 && ~ist0UTC
    error('请指定初始时刻（jd0或t0UTC）')
end
% t0TCB = aux.t0UTC;
% jd0 = juliandate(t0TCB(1), t0TCB(2), t0TCB(3), t0TCB(4), t0TCB(5), t0TCB(6));%ICRF坐标系,day
% aux.jd0 = jd0;% day


end