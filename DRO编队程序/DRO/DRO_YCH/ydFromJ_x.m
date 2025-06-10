function yd = ydFromJ_x(x, J,mu)

yd = sqrt(x^2+(2*(1-mu))/sqrt((x+mu)^2)+2*mu/sqrt((x-1+mu)^2) - J);