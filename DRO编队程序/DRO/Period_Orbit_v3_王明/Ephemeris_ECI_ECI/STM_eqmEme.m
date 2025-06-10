function xxdot = STM_eqmEme(t , xx , aux)

xxdot = zeros(42,1);
if size(xx, 2)~=1
   xx = xx';
end
mu_earth = aux.planet.mu(3);
mu_moon = aux.planet.mu(10);
mu_sun = aux.planet.mu(11);

% s -> day
jdate = t / 86400;

% 卫星位置矢量
r_sc = xx(1:3);

% 卫星位置矢量模
rmag_sc = sqrt( xx(1)^2 + xx(2)^2 + xx(3)^2 );

% 太阳位置矢量
% rv_sun = ephEme(jdate, 11, 3 , aux.DE430);
rv_sun = ephEme_mex(jdate, 11, 3 , aux.DE430);
r_sun = rv_sun(1:3);
rmag_sun = sqrt( r_sun(1)^2 + r_sun(2)^2 + r_sun(3)^2 );

% 月球位置矢量
% rv_moon = ephEme(jdate, 10, 3 , aux.DE430);
rv_moon = ephEme_mex(jdate, 10, 3 , aux.DE430);
r_moon = rv_moon(1:3);
rmag_moon = sqrt( r_moon(1)^2 + r_moon(2)^2 + r_moon(3)^2 );

% 太阳指向卫星位置矢量
r_sun2sc = r_sc - r_sun;
rmag_sun2sc = sqrt( r_sun2sc(1)^2 + r_sun2sc(2)^2 + r_sun2sc(3)^2 );

% 卫星指向月球位置矢量
r_moon2sc = r_sc - r_moon;
rmag_moon2sc = sqrt( r_moon2sc(1)^2 + r_moon2sc(2)^2 + r_moon2sc(3)^2 );

% --------------------- 计算加速度 ---------------------
% 加速度地球部分
acc_earth = - mu_earth * r_sc / rmag_sc^3;

% 加速度月球部分
acc_moon = - mu_moon * (r_moon / rmag_moon^3 + r_moon2sc / rmag_moon2sc^3);

% 加速度太阳部分
acc_sun = - mu_sun * (r_sun / rmag_sun^3 + r_sun2sc / rmag_sun2sc^3);

% 行星加速度
acc_planet = acc_earth + acc_moon + acc_sun;

% 卫星加速度
acc_sc = [0 ; 0 ; 0];

% compute integration vector
xxdot(1:6) = [ xx(4)
    xx(5)
    xx(6)
    acc_planet(1) + acc_sc(1);
    acc_planet(2) + acc_sc(2);
    acc_planet(3) + acc_sc(3)];

A1=zeros(3);
A2=eye(3);
A4=zeros(3);
G = zeros(3);
Main_body =[3, 10, 11];
rp2sc = zeros(3,11);
rp2sc(:,3) = r_sc;
rp2sc(:,10) = r_moon2sc;
rp2sc(:,11) = r_sun2sc;
for i =1:length(Main_body)
        index = Main_body(i);
        G(1,1) = G(1,1)+(-aux.planet.mu(index)*norm(rp2sc(:,index))^2+3*aux.planet.mu(index)*rp2sc(1,index)^2)/norm(rp2sc(:,index))^5;
        G(2,2) = G(2,2)+(-aux.planet.mu(index)*norm(rp2sc(:,index))^2+3*aux.planet.mu(index)*rp2sc(2,index)^2)/norm(rp2sc(:,index))^5;
        G(3,3) = G(3,3)+(-aux.planet.mu(index)*norm(rp2sc(:,index))^2+3*aux.planet.mu(index)*rp2sc(3,index)^2)/norm(rp2sc(:,index))^5;
        G(1,2) = G(1,2)+3*aux.planet.mu(index)*rp2sc(1,index)*rp2sc(2,index)/norm(rp2sc(:,index))^5;
        G(1,3) = G(1,3)+3*aux.planet.mu(index)*rp2sc(1,index)*rp2sc(3,index)/norm(rp2sc(:,index))^5;
        G(2,3) = G(2,3)+3*aux.planet.mu(index)*rp2sc(2,index)*rp2sc(3,index)/norm(rp2sc(:,index))^5;
end
G(2,1) = G(1,2);
G(3,1) = G(1,3);
G(3,2) = G(2,3);
A = [A1,A2;G A4];
xxdot(7:42) = reshape(A*reshape(xx(7:42),6,6),36,1);
end

