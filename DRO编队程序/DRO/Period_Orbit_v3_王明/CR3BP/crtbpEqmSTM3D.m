
function dxx = crtbpEqmSTM3D(t , xx , aux)
% 圆型限制性三体问题状态方程和状态转移矩阵的42维动力学方程
%
% 例子：
% aux = crtbpAuxEM;
% crtbpMarkEM
% xx0 = [0.75 , 0 , 0.2 , 0.15 , 0 , 0]';
% yy0 = [xx0 ; reshape(eye(6) , 36 , 1)];
% options = odeset('Reltol' , 1e-13, 'AbsTol' , 1e-13);
% [tt , yy] = ode113(@crtbpEqmSTM3D , [0 , 1e2] , yy0 , options , aux.mu);
% plot3D(yy(: , 1 : 3) , 'b' , 1);
% 
% 2015/4/10
% Copyright(C) Chen Zhang
% -------------------------------------------------------------------------

mu = aux.mu;

% State vector at t
x = xx(1);
y = xx(2);
z = xx(3);
vx = xx(4);
vy = xx(5);
vz = xx(6);

% STM at t
Phi = reshape(xx(7 : 42) , 6 , 6);

r1cube = ((x + mu)^2 + y^2 + z^2)^(3/2);
r2cube = ((x + mu - 1)^2 + y^2 + z^2)^(3/2);
r1five = r1cube^(5/3);
r2five = r2cube^(5/3);

G11 = 1 - (1 - mu) / r1cube + 3 * (1 - mu) * (x + mu)^2 / r1five - mu / r2cube + 3 * mu * (x + mu - 1)^2 / r2five;
G12 = 3 * (1 - mu) * (x + mu) * y / r1five + 3 * mu * (x + mu - 1) * y / r2five;
G13 = 3 * (1 - mu) * (x + mu) * z / r1five + 3 * mu * (x + mu - 1) * z / r2five;
G22 = 1 - (1 - mu) / r1cube + 3 * (1 - mu) * y^2 / r1five - mu / r2cube + 3 * mu * y^2 / r2five;
G23 = 3 * (1 - mu) * y * z / r1five + 3 * mu * y * z / r2five;
G33 = -(1 - mu) / r1cube + 3 * (1 - mu) * z^2 / r1five - mu / r2cube + 3 * mu * z^2 / r2five;

% jacobian of dynamics
Df = [0 , 0 , 0 , 1 , 0 , 0;
    0 , 0 , 0 , 0 , 1 , 0;
    0 , 0 , 0 , 0 , 0 , 1;
    G11 , G12 , G13 , 0 , 2 , 0;
    G12 , G22 , G23 , -2 , 0 , 0;
    G13 , G23 , G33 , 0 , 0 , 0];
dPhi = Df * Phi;

dxx = zeros(42 , 1);
dxx(1 : 6) = [vx;
    vy;
    vz;
    x - (1 - mu) * (x + mu) / r1cube - mu * (x + mu - 1) / r2cube + 2 * vy;
    y - (1 - mu) * y / r1cube - mu * y / r2cube - 2 * vx;
    -(1 - mu) * z / r1cube - mu * z / r2cube];
dxx(7 : 42) = dPhi(:);

end
