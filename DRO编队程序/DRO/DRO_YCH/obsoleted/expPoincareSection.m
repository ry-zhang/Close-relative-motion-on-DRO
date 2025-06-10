% clear; clc;
% Poincare Section (PS) with the same engergy
% start from x axis, [x 0 0 yd]
%%
global mu
mu = 0.01215; % Earth-Moon system

%%
[L] = lagrange_points(mu);
r1 = sqrt((L(:,1)+mu).^2+L(:,2).^2);
r2 = sqrt((L(:,1)-1+mu).^2+L(:,2).^2);
JLp=  (L(:,1).^2+L(:,2).^2)+2*((1-mu)./r1+mu./r2); % Energy of the L-points

%%
J = 2.88; % Energy
DT = 300; % intergration time
xMat = -1.5:0.03:1.5;
%  xMat = 0:0.02:1;
PS = cell(length(xMat), 1); % storing the PS data
Pd = cell(length(xMat), 1); % storing the periods
maxDimPd = 0;
for ii = 1:length(xMat)
    x = xMat(ii);
    [x,ii]
    % chk - check if the initial value lies inside the zero-velocity curve
    chk =  x^2+(2*(1-mu))/sqrt((x+mu)^2)+2*mu/sqrt((x-1+mu)^2)-J;
    temp = []; % [x y xd yd]
    T_temp = [];
    
    if chk>0 % if good
        yd = ydInit(x, J); % caculate the initial yd
        X0 = [x, 0, 0, yd]; % ics
        opts = odeset('RelTol',1e-9,'AbsTol',1e-9, 'Events',@secYeq0);
        [T,X, Te, Xe, Ie] = ode113(@(t,x)pcr3bp(t,x),[0 DT], X0, opts);
        temp = [Xe Te];
        if(length(Te)>1)
            T_temp = Te(2:end)-Te(1:end-1);
            maxDimPd = max(maxDimPd, length(T_temp));
        end
    end
    PS{ii,1} = temp;
    Pd{ii,1} = T_temp;
    
end

pdMat = [];
for ii = 1: length(xMat)
    pdMat(:,ii) = [Pd{ii,1}; nan*ones(maxDimPd-length(Pd{ii,1}),1)];
end

%%
close all

figure(1)
for ii = 1: length(xMat)
    temp = PS{ii,1};
    if ~isempty(temp)
        plot(temp(:,1), temp(:,3),'.') ; hold on
    end
end
plot(L(1:3,1),zeros(3,1),'rs');
grid on
xlabel('x'); ylabel('xd');
xlim([-1 1.5]);
ylim([-10, 10]);


% figure(2)
% boxplot(pdMat,xMat)