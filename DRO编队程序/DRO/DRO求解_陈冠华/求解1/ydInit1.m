function yd = ydInit1(x, xd, aux)

mu = aux.mu;
J = aux.J;

y = 0;  
r1 = sqrt((x+mu)^2+y^2);
r2 = sqrt((x-1+mu)^2+y^2);
yd = sqrt(x^2+y^2+(2*(1-mu))/r1+2*mu/r2-xd^2 - J);

