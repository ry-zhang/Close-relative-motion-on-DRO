function [rvRot,aNew]  = T_ECJ2k2Rot(jd_Mtx, rvJ2k_Mtx, a_j2k_Mtx, CoeffMat, targetFrame)
%%  将轨迹（Mtx）从地心惯性系ECJ2k-->其他旋转坐标系
%  targetFrame - 目标坐标系名称
% 大天体都在左边
% 输入：
% jd_Mtx       -  儒略日（day）
% rfJ2k_Mtx    -  轨迹rv
% CoeffMat     -  星历
% targetFrame  -  坐标系（字符串格式）
row = length(jd_Mtx);
rvRot = zeros(row,6);
aNew = zeros(row,3);
for ii = 1:row
    [rvRot(ii,:),aNew(ii,:)] = ECJ2k2Rot(jd_Mtx(ii), rvJ2k_Mtx(ii,:), a_j2k_Mtx(ii,:)', CoeffMat, targetFrame);
end
end

function [rvRot,aNew]  = ECJ2k2Rot(jd, rvJ2k, aJ2k, CoeffMat, targetFrame)
%%  地心惯性系ECJ2k-->其他坐标系
%  targetFrame - 目标坐标系名称
% 大天体都在左边

%  DE430 structure
%        1 = mercury           8 = neptune
%        2 = venus              9 = pluto
%        3 = EM                10 = moon(geo)
%        4 = mars             11 = sun
%        5 = jupiter           12 = nutation
%        6 = saturn           13 = moon lib
%        7 = uranus

rvJ2k = reshape(rvJ2k, 3, 2);
djd = 1e-4;

pEarthRV  = JPL_Eph_DE430_PosVel(jd,  3, CoeffMat); % 地球
pMoonRV  = JPL_Eph_DE430_PosVel(jd,  10, CoeffMat); % 月球
pSunRV  = JPL_Eph_DE430_PosVel(jd,  11, CoeffMat); % 太阳

pEarthRV_next  = JPL_Eph_DE430_PosVel(jd+djd,  3, CoeffMat); % 地球
pMoonRV_next  = JPL_Eph_DE430_PosVel(jd+djd,  10, CoeffMat); % 月球
pSunRV_next  = JPL_Eph_DE430_PosVel(jd+djd,  11, CoeffMat); % 太阳
pEarthRV_pre  = JPL_Eph_DE430_PosVel(jd-djd,  3, CoeffMat); % 地球
pMoonRV_pre  = JPL_Eph_DE430_PosVel(jd-djd,  10, CoeffMat); % 月球
pSunRV_pre  = JPL_Eph_DE430_PosVel(jd-djd,  11, CoeffMat); % 太阳

aEarth = (pEarthRV_next(:, 2)-pEarthRV_pre(:, 2))/djd/2; 
aMoon = (pMoonRV_next(:, 2)-pMoonRV_pre(:, 2))/djd/2;
aSun = (pSunRV_next(:, 2)-pSunRV_pre(:, 2))/djd/2;

adotEarth = (pEarthRV_next(:, 2)+pEarthRV_pre(:, 2)-2*pEarthRV(:, 2))/djd^2;
adotMoon = (pMoonRV_next(:, 2)+pMoonRV_pre(:, 2)-2*pMoonRV(:, 2))/djd^2;
adotSun = (pSunRV_next(:, 2)+pSunRV_pre(:, 2)-2*pSunRV(:, 2))/djd^2;

if strcmp(targetFrame, 'ECEMR')   %地心地月旋转系
    p1RV = pEarthRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pEarthRV; % 原点天体
    a1 = aEarth;  a2 = aMoon;  ao = aEarth;
    adot1 = adotEarth; adot2 = adotMoon;
elseif strcmp(targetFrame, 'MCEMR')  %月心地月旋转系
    p1RV = pEarthRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pMoonRV; % 原点天体
    a1 = aEarth; a2 = aMoon; ao = aMoon;
    adot1 = adotEarth; adot2 = adotMoon;
elseif strcmp(targetFrame, 'ECSER')  %地心日地旋转系
    p1RV = pSunRV; % 天体1
    p2RV = pEarthRV; % 天体2
    poRV = pEarthRV; % 原点天体
    a1 = aSun; a2 = aEarth; ao = aEarth;
    adot1 = adotSun; adot2 = adotEarth;
elseif strcmp(targetFrame, 'MCSMR')  %月心日月旋转系
    p1RV = pSunRV; % 天体1
    p2RV = pMoonRV; % 天体2
    poRV = pMoonRV; % 原点天体
    a1 = aSun; a2 = aMoon; ao = aMoon;
    adot1 = adotSun; adot2 = adotMoon;
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
rn = norm(r);
% vn = norm(v);
hn = norm(h);

iUnit = r/rn;
kUnit = h/hn;
jUnit = cross(kUnit, iUnit);
rotMat = [iUnit, jUnit, kUnit]';

a = (a2 - a1)/86400^2;
w = h/rn^2 + dot(a,h)/hn^2*r;
% w = h/rn^2;

adot = (adot2 - adot1)/86400^3;
hdot = cross(r,a);
hndot = dot(hdot,kUnit);
rndot = dot(r,v)/rn;
wdot = - 2*rndot/rn^3*h + hdot/rn^2 + ...
    - 2*hndot/hn^3*dot(a,h)*r + (dot(adot,h)+dot(a,hdot))/hn^2*r + dot(a,h)/hn^2*v;


dr_I = rvJ2k(:,1) - (poRV(:,1)-pEarthRV(:,1));
dv_I = rvJ2k(:,2) - (poRV(:,2)-pEarthRV(:,2))/86400;
    
da_I = aJ2k - (ao-aEarth)/86400^2;

rvNew(:, 1) =  rotMat*(dr_I);
ve = cross(w, dr_I);
rvNew(:, 2) =  rotMat*(dv_I-ve);
rvRot = rvNew(:);

aNew = rotMat*(da_I-cross(wdot,dr_I)-2*cross(w,dv_I)+cross(w,cross(w,dr_I)));

end