function rvJ2k = T_Rot2ECJ2k(jd_Mtx, rvRot_Mtx, CoeffMat, fromFrame)
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
row = length(jd_Mtx);
rvJ2k = zeros(row,6);
for ii = 1:row
    rvJ2k(ii,:) = Rot2ECJ2k(jd_Mtx(ii), rvRot_Mtx(ii,:), CoeffMat, fromFrame);
end

end

function  rvJ2k = Rot2ECJ2k(jd, rvRot, CoeffMat, fromFrame)
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

% rvRot - column vector

rvRot = reshape(rvRot, 3, 2);
djd = 1e-6;

pEarthRV  = JPL_Eph_DE430_PosVel(jd,  3, CoeffMat); % 地球
pMoonRV  = JPL_Eph_DE430_PosVel(jd,  10, CoeffMat); % 月球
pSunRV  = JPL_Eph_DE430_PosVel(jd,  11, CoeffMat); % 太阳

pEarthRV_next  = JPL_Eph_DE430_PosVel(jd+djd,  3, CoeffMat); % 地球
pMoonRV_next  = JPL_Eph_DE430_PosVel(jd+djd,  10, CoeffMat); % 月球
pSunRV_next  = JPL_Eph_DE430_PosVel(jd+djd,  11, CoeffMat); % 太阳

aEarth = (pEarthRV_next(:, 2)-pEarthRV(:, 2))/djd; 
aMoon = (pMoonRV_next(:, 2)-pMoonRV(:, 2))/djd;
aSun = (pSunRV_next(:, 2)-pSunRV(:, 2))/djd;

if strcmp(fromFrame, 'ECEMR')   %地心地月旋转系
    p1RV = pEarthRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pEarthRV; % 原点天体
    a1 = aEarth;  a2 = aMoon;
elseif strcmp(fromFrame, 'MCEMR')  %月心地月旋转系
    p1RV = pEarthRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pMoonRV; % 原点天体
    a1 = aEarth; a2 = aMoon;
elseif strcmp(fromFrame, 'ECSER')  %地心日地旋转系
    p1RV = pSunRV; % 天体1
    p2RV = pEarthRV; % 天体2
    poRV = pEarthRV; % 原点天体
    a1 = aSun; a2 = aEarth;
elseif strcmp(fromFrame, 'MCSMR')  %月心日月旋转系
    p1RV = pSunRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pMoonRV; % 原点天体
    a1 = aSun; a2 = aMoon;
else
    error('Wrong target frame');
end


rp1 =  p1RV(:, 1);
vp1 =  p1RV(:, 2);
rp2 =  p2RV(:, 1);
vp2 =  p2RV(:, 2);

r = rp2 - rp1;
v = (vp2 - vp1)/86400; % 注意DE430单位, km/day
h = cross(r, v);
a = (a2 - a1)/86400^2;
w = h/norm(r)^2 + dot(a,h)/norm(h)^2*r;
% w = h/norm(r)^2;

iUnit = r/norm(r);
kUnit = h/norm(h);
jUnit =  cross(kUnit, iUnit);
rotMat = [iUnit, jUnit, kUnit];
 
% rot origin wrp to Earth
rOE = poRV(:,1)-pEarthRV(:,1);
vOE = (poRV(:,2)-pEarthRV(:,2))/86400;

rvNew(:, 1)  = rotMat*rvRot(:, 1)+ (poRV(:,1)-pEarthRV(:,1));

ve = cross(w,  rvNew(:, 1) -rOE);
rvNew(:, 2)  = rotMat*rvRot(:, 2)+ve+vOE;

rvJ2k = rvNew(:);
end