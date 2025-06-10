function [Df, dxx] = phiPcr3bpAug(t,  xx , auxdata)

dim = 4;
mu = auxdata.mu;

x = xx(1); y = xx(2); vx = xx(3); vy = xx(4);
phi = reshape(xx(dim+1:end), dim, dim);

mu1 = 1-mu; % mass of larger  primary (nearest origin on left)
mu2 =   mu; % mass of smaller primary (furthest from origin on right)

r2= (x+mu2)^2 + y^2;      % r: distance to m1, LARGER MASS
R2= (x-mu1)^2 + y^2;      % R: distance to m2, smaller mass

r3= r2^1.5;
r5= r2^2.5;
R3= R2^1.5;
R5= R2^2.5;

Ux = - x + mu1*(x+mu2)/r3 + mu2*(x-mu1)/R3 ;
Uy = - y + mu1* y     /r3 + mu2* y     /R3 ;

Uxx = -1+(mu1/r3)*(1-(3*(x+mu2)^2/r2))+(mu2/R3)*(1-(3*(x-mu1)^2/R2)) ;
Uyy = -1+(mu1/r3)*(1-(3* y     ^2/r2))+(mu2/R3)*(1-(3* y     ^2/R2)) ;
Uxy =   -(mu1/r5)*    3* y*(x+mu2) -(mu2/R5)*    3* y*(x-mu1)  ;


Df = [  0     0    1    0 ;
          0     0    0    1 ;
	-Uxx  -Uxy   0    2 ;
        -Uxy  -Uyy  -2    0 ];


phiDot = Df*phi;

dxx = zeros(dim+dim^2,1);

dxx(1) = vx;
dxx(2) = vy;
dxx(3) = 2*vy - Ux ;
dxx(4) =-2*vx - Uy ;

dxx(dim+1:end) = phiDot(:);

end