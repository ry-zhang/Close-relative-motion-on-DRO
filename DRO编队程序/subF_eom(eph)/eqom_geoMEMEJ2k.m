function ydot = eqom_geoMEMEJ2k(t, xx, aux)
% --------------------- 功能--------------
% 行星际飞行运动积分
% 坐标系：geo MEME J2k，地心J2K平赤道平春分点坐标系
%  主天体：地球
%  主要摄动：太阳、月球
%  摄动：行星1-2，4-9/3可选, 3=EM
% --------------------------输入-------------
%  t = current simulation time (day)
% output
%  ydot = first order equations of motion


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  transpose xx to a coloum vector
if size(xx, 2)~=1
    xx = xx';
end

% extract data
jd0 = aux.jd0;  %参考儒略日
xmu = aux.xmu;
C_Mat = aux.C_Mat;
% ICRF_2_MeanEclpJ2k = aux.ICRF_2_MeanEclpJ2k;
ICRF_2_J2K =  aux.ICRF_2_J2K;
threeBody = aux.threeBodyList;
if any(threeBody>=10) || any(threeBody==3)
    error('3 = EM, 月球已经考虑过了')
end

mainBody =[3, 10, 11];

% current julian date (relative to tcm event)
jdate = jd0 + t / 86400;

% distance from Earth to spacecraft
rEarth2sc = norm(xx(1:3));
rrEarth2sc = -xmu(3) / rEarth2sc^3; % mu/r^3

pMat = zeros(3,11);
for ii = 1:11 %所有天体的ssb ICRF位置
    if any(ii == mainBody) || any(ii == threeBody)
        pv  =  JPL_Eph_DE430_PosVel(jdate,  ii,  C_Mat);
    else
        pv = zeros(3,1);
    end
    pMat(:, ii)= pv(:,1);
end

rp = zeros(3,11);
rp2sc = zeros(3,11);
for ii = 1:11 %所有天体相对earth和sc的位置， TEME系
    rp(:, ii) = ICRF_2_J2K * (pMat(:, ii)-pMat(:,3)); %天体相对于earth
    rp2sc(:, ii)  = xx(1:3) - rp(:, ii); %天体相对于sc（航天器相对于天体吧）
end

accp = zeros(3,1);
for ii = 10:11%计算主要摄动：月球 太阳
    d = rp2sc(:, ii) ;%月球/太阳→航天器
    rho = rp(:, ii);%地球→月球/太阳
    accp = accp - xmu(ii) * (d*norm(d)^(-3)+rho*norm(rho)^(-3));
end

%↓（计算行星摄动）
for ii = 1:length(threeBody)
    idPlanet = threeBody(ii);
    d = rp2sc(:, idPlanet) ;
    rho = rp(:, idPlanet);
    accp = accp - xmu(idPlanet) * (d*norm(d)^(-3)+rho*norm(rho)^(-3));
end

% compute integration vector
ydot = [ xx(4:6)
    accp  + xx(1:3)  * rrEarth2sc];
end

