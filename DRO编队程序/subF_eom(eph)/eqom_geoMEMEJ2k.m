function ydot = eqom_geoMEMEJ2k(t, xx, aux)
% --------------------- ����--------------
% ���Ǽʷ����˶�����
% ����ϵ��geo MEME J2k������J2Kƽ���ƽ���ֵ�����ϵ
%  �����壺����
%  ��Ҫ�㶯��̫��������
%  �㶯������1-2��4-9/3��ѡ, 3=EM
% --------------------------����-------------
%  t = current simulation time (day)
% output
%  ydot = first order equations of motion


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  transpose xx to a coloum vector
if size(xx, 2)~=1
    xx = xx';
end

% extract data
jd0 = aux.jd0;  %�ο�������
xmu = aux.xmu;
C_Mat = aux.C_Mat;
% ICRF_2_MeanEclpJ2k = aux.ICRF_2_MeanEclpJ2k;
ICRF_2_J2K =  aux.ICRF_2_J2K;
threeBody = aux.threeBodyList;
if any(threeBody>=10) || any(threeBody==3)
    error('3 = EM, �����Ѿ����ǹ���')
end

mainBody =[3, 10, 11];

% current julian date (relative to tcm event)
jdate = jd0 + t / 86400;

% distance from Earth to spacecraft
rEarth2sc = norm(xx(1:3));
rrEarth2sc = -xmu(3) / rEarth2sc^3; % mu/r^3

pMat = zeros(3,11);
for ii = 1:11 %���������ssb ICRFλ��
    if any(ii == mainBody) || any(ii == threeBody)
        pv  =  JPL_Eph_DE430_PosVel(jdate,  ii,  C_Mat);
    else
        pv = zeros(3,1);
    end
    pMat(:, ii)= pv(:,1);
end

rp = zeros(3,11);
rp2sc = zeros(3,11);
for ii = 1:11 %�����������earth��sc��λ�ã� TEMEϵ
    rp(:, ii) = ICRF_2_J2K * (pMat(:, ii)-pMat(:,3)); %���������earth
    rp2sc(:, ii)  = xx(1:3) - rp(:, ii); %���������sc�����������������ɣ�
end

accp = zeros(3,1);
for ii = 10:11%������Ҫ�㶯������ ̫��
    d = rp2sc(:, ii) ;%����/̫����������
    rho = rp(:, ii);%���������/̫��
    accp = accp - xmu(ii) * (d*norm(d)^(-3)+rho*norm(rho)^(-3));
end

%�������������㶯��
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

