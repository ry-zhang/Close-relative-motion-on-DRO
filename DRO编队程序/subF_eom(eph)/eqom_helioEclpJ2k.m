function ydot = eqom_helioEclpJ2k(t, xx, aux)
% --------------------- 功能--------------
% 行星际飞行运动积分
% 坐标系：sun EclpJ2k，日心平黄道惯性系
%  主天体：太阳
%  摄动：行星1-9可选, 3=EM
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
ICRF_2_MeanEclpJ2k = aux.ICRF_2_MeanEclpJ2k;
threeBody = aux.threeBodyList;
if any(threeBody>=10)
    error('3 = EM, 月球已经考虑过了')
end


% current julian date (relative to tcm event)
jdate = jd0 + t / 86400;

% distance from sun to spacecraft
rsun2sc = norm(xx(1:3));
rrsun2sc = -xmu(11) / rsun2sc^3; % mu/r^3


for ii = 1:11 %所有天体的ssb ICRF位置
    pv  =  JPL_Eph_DE430_PosVel(jdate,  ii,  C_Mat);
    pMat(:, ii)= pv(:,1);
end

for ii = 1:11 %所有天体相对sun和sc的位置， Eclp系
    rp(:, ii) = ICRF_2_MeanEclpJ2k * (pMat(:, ii)-pMat(:,11)); %天体相对于太阳
    rp2sc(:, ii)  = xx(1:3) - rp(:, ii);
end


% compute planetary perturbations
accp = zeros(3,1);
for ii = 1:length(threeBody)
    idPlanet = threeBody(ii);
    d = rp2sc(:, idPlanet) ;
    rho = rp(:, idPlanet);
    accp = accp - xmu(idPlanet) * (d*norm(d)^(-3)+rho*norm(rho)^(-3));
end

% compute integration vector
ydot = [ xx(4:6)
    accp  + xx(1:3)  * rrsun2sc];
end


