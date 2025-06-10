function yd = ydInit(x, J)

global mu
y = 0; xd = 0;
r1 = sqrt((x+mu)^2+y^2);
r2 = sqrt((x-1+mu)^2+y^2);
yd = sqrt(x^2+y^2+(2*(1-mu))/r1+2*mu/r2-xd^2 - J);