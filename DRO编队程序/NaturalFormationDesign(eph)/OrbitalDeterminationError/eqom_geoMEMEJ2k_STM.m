function dxdx = eqom_geoMEMEJ2k_STM(t,xx,aux)
% --------------------- 功能--------------
% 行星际飞行运动积分
% 坐标系：geo MEME J2k，地心J2K平赤道平春分点坐标系
%  主天体：地球
%  主要摄动：太阳、月球
%  摄动：行星1-2，4-9/3可选, 3=EM
% --------------------------输入-------------
%  t = current simulation time (day)
% output
%  dxdx = first order equations of motion


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  transpose xx to a coloum vector


dxdx = zeros(42,1);
if size(xx, 2)~=1
    xx = xx';
end


jd0 = aux.jd0;  % 参�?儒略�?
xmu = aux.xmu;
C_Mat = aux.C_Mat;
ICRF_2_MeanEclpJ2k = aux.ICRF_2_MeanEclpJ2k;
ICRF_2_J2K =  aux.ICRF_2_J2K;
threeBody = aux.threeBodyList;
if any(threeBody>=10) || any(threeBody==3)
   error('3 = EM, 月球已经考虑过了')
end
jdate = jd0 + t / 86400;
rEarth2sc = norm(xx(1:3));  % 航天器在地心J2k坐标系中的位�?
rrEarth2sc = -xmu(3) / rEarth2sc^3; % mu/r^3
for ii = 1:11 %�?��天体的ssb ICRF位置
    pv  =  JPL_Eph_DE430_PosVel(jdate,  ii,  C_Mat);  %%这个函数输出位置和�?�?
    pMat(:, ii)= pv(:,1);
end
for ii = 1:11 %�?��天体相对earth和sc的位置， TEME�?
    rp(:, ii) = ICRF_2_J2K * (pMat(:, ii)-pMat(:,3)); %天体相对于earth
    rp2sc(:, ii)  = xx(1:3) - rp(:, ii); %天体相对于sc
end

accp = zeros(3,1);
for ii = 10:11
    d = rp2sc(:, ii) ;
    rho = rp(:, ii);
    accp = accp - xmu(ii) * (d*norm(d)^(-3)+rho*norm(rho)^(-3));
end

for ii = 1:length(threeBody)
    idPlanet = threeBody(ii);
    d = rp2sc(:, idPlanet) ;
    rho = rp(:, idPlanet);
    accp = accp - xmu(idPlanet) * (d*norm(d)^(-3)+rho*norm(rho)^(-3));
end

% compute integration vector
dxdx(1:6) = [ xx(4:6);
             accp  + xx(1:3)  * rrEarth2sc];
 A1=zeros(3);
 A2=eye(3);
 A4=zeros(3);
 G = zeros(3);
 Main_body =[3, 10, 11];

for i =1:length(Main_body)
    index = Main_body(i);
    G(1,1) = G(1,1)+(-xmu(index)*norm(rp2sc(:,index))^2+3*xmu(index)*rp2sc(1,index)^2)/norm(rp2sc(:,index))^5;
    G(2,2) = G(2,2)+(-xmu(index)*norm(rp2sc(:,index))^2+3*xmu(index)*rp2sc(2,index)^2)/norm(rp2sc(:,index))^5;
    G(3,3) = G(3,3)+(-xmu(index)*norm(rp2sc(:,index))^2+3*xmu(index)*rp2sc(3,index)^2)/norm(rp2sc(:,index))^5;
    G(1,2) = G(1,2)+3*xmu(index)*rp2sc(1,index)*rp2sc(2,index)/norm(rp2sc(:,index))^5;
    G(1,3) = G(1,3)+3*xmu(index)*rp2sc(1,index)*rp2sc(3,index)/norm(rp2sc(:,index))^5;
    G(2,3) = G(2,3)+3*xmu(index)*rp2sc(2,index)*rp2sc(3,index)/norm(rp2sc(:,index))^5;
end
G(2,1) = G(1,2);
G(3,1) = G(1,3);
G(3,2) = G(2,3);
A = [A1,A2;G A4];
dxdx(7:42) = reshape(A*reshape(xx(7:42),6,6),36,1);

end
