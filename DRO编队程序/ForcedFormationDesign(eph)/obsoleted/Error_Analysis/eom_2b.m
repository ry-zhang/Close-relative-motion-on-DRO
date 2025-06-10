function rrdot = eom_2b(t,y)
mu = 398600; % km3/s2
rrdot = zeros(6,1);
% position and velocity
rx = y(1);
ry = y(2);
rz = y(3);
vx = y(4);
vy = y(5);
vz = y(6);
% differnetial equations
r = sqrt(rx^2 + ry^2 + rz^2);
r3 = r^3;
rrdot(1) = vx;
rrdot(2) = vy;
rrdot(3) = vz;
rrdot(4) = mu*(-rx/r3);
rrdot(5) = mu*(-ry/r3);
rrdot(6) = mu*(-rz/r3);
