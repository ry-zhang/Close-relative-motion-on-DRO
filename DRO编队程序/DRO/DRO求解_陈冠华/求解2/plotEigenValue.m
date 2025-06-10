% 画画专用
clear; clc

mu = 0.01215;
aux.mu = mu;
aux.dim = 4; %平面问题;

load('DRO_1po_2.4_3.1_0.005_.mat')
% J_x1(:,[1 5]) = dro3p(:,[5 1]);


phi0 = reshape(eye(4),1,16);
opts0 = odeset('RelTol',1e-13,'AbsTol',1e-13, 'Events',@secYeq0);
DT = 20;

yd = []; xx = [];
phiDot = cell(length(dro3p),1);
lambda = [];lam1=[]; lam = [];

for ii = 1:length(dro3p)
    aux.J = dro3p(ii,5);
    yd(ii) = ydInit2(dro3p(ii,1), 0, aux);
    xx(ii,:) = [dro3p(ii,1), 0, 0, yd(ii), phi0];
[~,~, TeVar, XeVar, ~] = ode113(@(t,x)odePcr3bpAug(t ,x, aux),[0 DT], xx(ii,:), opts0);
    phiDot{ii} = reshape(XeVar(end,5:end),4,4);
    lambda(ii,:) = eig(phiDot{ii});
    lam(ii,:) = abs(lambda(ii,:));
    theta(ii,:) = angle(lambda(ii,:))*180/pi;
end
J_theta =[ dro3p(:,5), theta];
J_lambda = [ dro3p(:,5), lambda];
save('J_x0_λ&θ_.mat', 'J_theta','J_lambda');


figure(1)
plot(J_lambda(:,1), J_lambda(:,2),'.',...
    J_lambda(:,1), J_lambda(:,3),'.',...
    J_lambda(:,1), J_lambda(:,4),'.',...
    J_lambda(:,1), J_lambda(:,5),'.')
hold on; grid on
xlabel('J');ylabel('λ')
xlim([2.4,3.1])
title('1-p')

figure(2)
plot(J_theta(:,1), J_theta(:,2),'.',...
    J_theta(:,1), J_theta(:,3),'.',...
    J_theta(:,1), J_theta(:,4),'.',...
    J_theta(:,1), J_theta(:,5),'.')
hold on; grid on
% legend('1-p')
xlabel('J');ylabel('θ')
xlim([2.4,3.1])
title('1-p')