function posVel = jplEph(JD,  eObj, CoeffMat)
% compute ssb pos and vel
% ----------input-----------
% JD -- Julian date
% eObj -- 1~9=Mercury to Pluto
% eObj -- 10=moon, 11 = sun
%
% ----------output --------
% % from ssb to celetrial body
% unit - km, day
% 3*2 mat
%--------------
% Copyright(C) 2019 by Hao Zhang @CSU,CAS
% hao.zhang.zhr@gmail.com
%


%% code starts here
% main
%  DE430 structure
%        1 = mercury           8 = neptune
%        2 = venus              9 = pluto
%        3 = EM                10 = moon(geo)
%        4 = mars             11 = sun
%        5 = jupiter           12 = nutation
%        6 = saturn           13 = moon lib
%        7 = uranus

header = [  1       2       3       4       5        6       7        8       9       10     11     12   13;            % legend
    3     171   231   309   342   366   387   405   423   441   753   819   899;         % start
    14    10    13      11      8        7        6       6        6      13     11     10     10;           % # of coeff
    4      2      2         1      1        1        1       1        1       8        2       4       4];            % # of subinterval

EMRAT = 81.30056907419062; % DE430  % m_Earth/m_Moon
EMRAT1 = 1/(1+EMRAT);

%% get the coeffs by checking  JD
ii = find(CoeffMat(: , 1) <= JD & JD < CoeffMat(: , 2) , 1 , 'first');
PCtemp = CoeffMat(ii,:); % the entire row
t1 = PCtemp(1); % JD at start of interval
dt = JD - t1;

%% get the coeffs of the planet
% [nc*3, nc*3, nc*3,...]

icStart = header(2, eObj);
nc = header(3, eObj);
ns =header(4, eObj);
interval_length = 32/ ns;
jIntv = fix(dt/interval_length)+1;
offset = dt-(jIntv-1)*interval_length;
tau = 2.0 * offset / interval_length - 1.0;

coeff = getCoeff(icStart, nc, ns, jIntv,PCtemp);
rv = cheb_Interp(coeff, tau, nc);
rv(:,2) = rv(:,2)*2 /interval_length;

if eObj == 3 % Earth
    icStart = header(2, 10);
    nc = header(3, 10);
    ns =header(4, 10);
    interval_length = 32/ ns;
    jIntv = fix(dt/interval_length)+1;
    offset = dt-(jIntv-1)*interval_length;
    tau = 2.0 * offset / interval_length - 1.0;
    
    coeff = getCoeff(icStart, nc, ns, jIntv,PCtemp);
    rvMnGeo = cheb_Interp(coeff, tau, nc); % Moon @geo
    rvMnGeo(:,2) = rvMnGeo(:,2)*2 /interval_length;
    
    % Earth from Solar Bary
    rv = rv-EMRAT1*rvMnGeo;
end

if eObj == 10 % Moon
    icStart = header(2, 3);
    nc = header(3, 3);
    ns =header(4, 3);
    interval_length = 32/ ns;
    jIntv = fix(dt/interval_length)+1;
    offset = dt-(jIntv-1)*interval_length;
    tau = 2.0 * offset / interval_length - 1.0;
    
    coeff = getCoeff(icStart, nc, ns, jIntv,PCtemp);
    rvEM = cheb_Interp(coeff, tau, nc); % EM @ssb
    rvEM(:,2) = rvEM(:,2)*2 /interval_length;
    
    rv = (1-EMRAT1)*rv+rvEM;
end

posVel = rv;
end
