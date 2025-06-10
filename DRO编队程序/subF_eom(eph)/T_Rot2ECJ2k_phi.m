function  [Phi_Rot2ECJ2k,Phidot_Rot2ECJ2k] = T_Rot2ECJ2k_phi(jd, CoeffMat, fromFrame)
%%  其他坐标系-->地心惯性系ECJ2k
% 'ECEMR'
% 'MCEMR'
% 'ECSER'
% 'MCSMR'

%  fromFrame - 原坐标系名称
% 大天体都在左边

%  DE430 structure
%        1 = mercury           8 = neptune
%        2 = venus              9 = pluto
%        3 = EM                10 = moon(geo)
%        4 = mars             11 = sun
%        5 = jupiter           12 = nutation
%        6 = saturn           13 = moon lib
%        7 = uranus


pEarthRV  = JPL_Eph_DE430_PosVelAccDacc(jd,  3, CoeffMat); % 地球
pMoonRV  = JPL_Eph_DE430_PosVelAccDacc(jd,  10, CoeffMat); % 月球
pSunRV  = JPL_Eph_DE430_PosVelAccDacc(jd,  11, CoeffMat); % 太阳

if strcmp(fromFrame, 'ECEMR')   %地心地月旋转系
    p1RV = pEarthRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pEarthRV; % 原点天体
elseif strcmp(fromFrame, 'MCEMR')  %月心地月旋转系
    p1RV = pEarthRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pMoonRV; % 原点天体
elseif strcmp(fromFrame, 'ECSER')  %地心日地旋转系
    p1RV = pSunRV; % 天体1
    p2RV = pEarthRV; % 天体2
    poRV = pEarthRV; % 原点天体
elseif strcmp(fromFrame, 'MCSMR')  %月心日月旋转系
    p1RV = pSunRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pMoonRV; % 原点天体
else
    error('Wrong target frame');
end

rp1 =  p1RV(:, 1);
vp1 =  p1RV(:, 2);
ap1 = p1RV(:, 3);
dap1 = p1RV(:, 4);
rp2 =  p2RV(:, 1);
vp2 =  p2RV(:, 2);
ap2 = p2RV(:, 3);
dap2 = p2RV(:, 4);

r = rp2 - rp1;
v = (vp2 - vp1)/86400; % 注意DE430单位, km/day
h = cross(r, v);
rn = norm(r);
% vn = norm(v);
hn = norm(h);

a = (ap2 - ap1)/86400^2;
adot = (dap2 - dap1)/86400^3;
% w = h/rn^2 + dot(a,h)/hn^2*r;
% w = h/norm(r)^2; 

kUnit = h/hn;
iUnit = r/rn;
jUnit = cross(kUnit,iUnit);
rotMat = [iUnit, jUnit, kUnit];

iUnit_dot = dot(v,jUnit) / rn * jUnit;
kUnit_dot = - dot(a,kUnit) * rn / hn * jUnit;
jUnit_dot = rn/hn * dot(a,kUnit)*kUnit - 1/rn * dot(v,jUnit)*iUnit;
rotMat_dot = [iUnit_dot, jUnit_dot, kUnit_dot];

Phi_Rot2ECJ2k = [rotMat, zeros(3); rotMat_dot, rotMat];

if nargout == 2
    iUnit_ddot = dot(a,jUnit)*jUnit/rn + dot(v,jUnit_dot)*jUnit/rn - dot(r,v)/rn^2*jUnit_dot + dot(v,jUnit)/rn*jUnit_dot;
    kUnit_ddot = (1/hn^2-dot(r,v)/rn^2)*kUnit_dot - rn/hn*(dot(adot,kUnit)+dot(a,kUnit_dot))*jUnit - rn/hn*dot(a,kUnit)*jUnit_dot;
    jUnit_ddot = dot(kUnit_ddot,iUnit) + 2*dot(kUnit_dot,iUnit_dot) + dot(kUnit,iUnit_ddot);
    rotMat_ddot = [iUnit_ddot, jUnit_ddot, kUnit_ddot];

    Phidot_Rot2ECJ2k = [rotMat_dot, zeros(3); rotMat_ddot, rotMat_dot];
end
end