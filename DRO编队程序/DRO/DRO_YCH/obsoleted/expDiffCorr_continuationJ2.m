clear; clc;
% find periodic orbit

global mu
mu = 0.01215; % 20200531

L = lagrange_points(mu);
r1 = sqrt((L(:,1)+mu).^2+L(:,2).^2);
r2 = sqrt((L(:,1)-1+mu).^2+L(:,2).^2);
JLp=  (L(:,1).^2+L(:,2).^2)+2*((1-mu)./r1+mu./r2);

%%
% x0 = 0.83; 
% JJ = [2.7:0.025:3.3];

% x0 = 0.3; 
% JJ = [2.965];

% 1/2的轨道周期
JJ = linspace(2.93052374520892,2.93052374520893,5);
x0 = 0.83;
T_DRO_desired = 1/2;

% pi/8的轨道周期
% JJ = linspace(2.95394755678542,2.95394755678543,5);
% x0 = 0.83;
% T_DRO_desired = pi/8;

% 1/3的轨道周期
% JJ = linspace( 2.97001118561685,2.97001118561688,5); %
% x0 = 0.85;
% T_DRO_desired = 1/3;

% 1/4的轨道周期
% JJ = linspace(3.00048914754381,3.00048914754383,5);
% x0 = 0.885;
% T_DRO_desired = 1/4;



plot_flag = 0;
state_ini_all = zeros(length(JJ),4);
J_period_all = zeros(length(JJ),2);

%% 
xiter = x0;
Pd=[];
for jj = 1:length(JJ)
    J = JJ(jj);
    iter= 1;

    [y, tp, D] = mapErr(xiter, J);
    while norm(y)>1e-13 && iter <50
        xiter = (xiter- y/D);
        [y, tp, D] = mapErr(xiter, J);
        iter = iter +1;
    end
    xPdc = xiter;

    ydPdc = ydFromJ_x(xPdc, J);


    %%
    opts1 = odeset('RelTol',1e-13,'AbsTol',1e-20,'Events',@secYeq0);
    DT = 100;
    state_ini_all(jj,:) = [xPdc, 0, 0, ydPdc];
    sol = ode113(@(t,x)pcr3bp(t,x),[0 DT], [xPdc, 0, 0, ydPdc], opts1);
    T_temp = diff(sol.xe);
    delT = T_temp(4);

    [J,delT/(2*pi)-T_DRO_desired]
    J_period_all(jj,:) = [J,delT/(2*pi)];
    Pd = [Pd; delT];


%%
if plot_flag == 1
    figure(3)
    t_I = sol.x;
    X_I =  zeros(size(sol.y));
    for i = 1:length(t_I)
        t_epoch = t_I(i);
        X_I([1,2],i) = synodic2inertial(sol.y([1,2],i),t_epoch);
    end
    subplot(2,2,1)
    plot(sol.y(1,:),sol.y(2,:)); axis equal; hold on
    plot(-mu,0,'.b','MarkerSize',30); 
    plot(1-mu,0,'.g','MarkerSize',10,'Color',[150, 150, 150]/255); hold off
    title('synodic coordinate')
    set(gca,'FontSize',13)
%     xlim([-0.3 1,3])

    subplot(2,2,2)
    plot(sol.x,sol.y(1,:),sol.x,sol.y(2,:))
    grid on; legend('x','y')
    title('synodic coordinate')
    set(gca,'FontSize',13)
    
    subplot(2,2,3)
    plot(X_I(1,:),X_I(2,:)); axis equal
    title('inertial coordinate')
    set(gca,'FontSize',13)
    
    subplot(2,2,4)
    plot(sol.x,X_I(1,:),sol.x,X_I(2,:))
    grid on; legend('x','y')
    title('inertial coordinate')
    set(gca,'FontSize',13)
    
    pause(0.0001)
end
    
end
save DRO_all1 state_ini_all J_period_all
%%
% figure(4)
% Tfit = fit(JJ',Pd,'smoothingspline');
% plot(Tfit,JJ',Pd);
% xlabel('J')
% ylabel('delT')






