clear; clc;
% find periodic orbit

global mu
mu = 0.01211;

L = lagrange_points(mu);
r1 = sqrt((L(:,1)+mu).^2+L(:,2).^2);
r2 = sqrt((L(:,1)-1+mu).^2+L(:,2).^2);
JLp=  (L(:,1).^2+L(:,2).^2)+2*((1-mu)./r1+mu./r2);


plot(L(:,1), L(:,2), 's');hold on
grid on
axis equal
plot(-mu,0,'o');
plot(1-mu,0,'o');
xlim([-1.2,1.2]);
ylim([-1,1]);
xlabel('x/LU');
ylabel('y/LU');