clear; clc;
% trajectory
%%
global mu
mu = 0.01211;

%%
[L] = lagrange_points(mu);
r1 = sqrt((L(:,1)+mu).^2+L(:,2).^2);
r2 = sqrt((L(:,1)-1+mu).^2+L(:,2).^2);
JLp=  (L(:,1).^2+L(:,2).^2)+2*((1-mu)./r1+mu./r2);

%%
%X0 = [0.1, 0, 0, 3.968047179119731]; % row
%X0 = [-0.25, 0, 0, 2.281994924]; % period - 1
%X0 = [-0.4, 0,0, 1.446067607]; % period 9

% c = 2.9251
% X0 = [0.108, 0, 0, 3.648142618]; % period 7
% X0 = [0.22, 0, 0.9, 2.145273432]; % period 3
% X0 = [0.944, 0, 0, .5739754672]; % moon orbit
% X0 = [0.8, 0, 0, 0.5260658834480009]; % DRO
% X0 = [0.335, 0, 0, 1.681549191]; % stable orbit around earth
% X0 = [0.25, 0, 1, 1.858306638];  % period-3 stable
% X0 = [0.84945, 0, 0, 0.9886090734e-1];  %

% %% J = C3
% X0 = [0.1528, 0, 1.394, 2.660130848]; % period-3 stable
% X0 = [0.334, 0, 0, 1.686223579]; %  number 8 stable
% X0 = [0.5902, 0, 0, .8220688736]; %  triangle stable
% X0 = [0.894, 0, 0, .4732097315]; %  DRO stable
% X0 = [0.87, 0, 0.044, .4474860769]; %  DRO stable


%% J = 2.7
% X0 = [0.4935, 0, 0, 1.224845846]; %  DRO stable
% X0 = [0.5228, 0, 0, 1.148505061]; %  4PDRO stable
% X0 = [0.85, 0, 0.16, 1.237868905]; %  QDRO stable

% %% J=3.18
% X0 = [0.944, 0, 0,  0.5739754672]; % DRO
% X0 = [0.964, 0, 0.46, 0.7587133418]; % QDRO
% X0 = [1.01, 0, 0, .9319723204];

DT = 300;

%%
J = 2.9251;
% x0 = 0.2725; % stable orbit around earth
 x0 = 0.6; % -1.808


yd0 = ydInit(x0, J);
X0 = [x0, 0, 0, yd0];

%% classical intergrator
%opts = odeset('RelTol',1e-11,'AbsTol',1e-11);
opts = odeset('RelTol',1e-12,'AbsTol',1e-12, 'Events',@secYeq0);

[T,X, ] = ode113(@(t,x)pcr3bp(t,x),[0 DT], X0, opts);

r1 = sqrt((X(:,1)+mu).^2+X(:,2).^2);
r2 = sqrt((X(:,1)-1+mu).^2+X(:,2).^2);

Jode =  (X(:,1).^2+X(:,2).^2)+2*((1-mu)./r1+mu./r2)-(X(:,3).^2+X(:,4).^2);

close all

figure(1)
plot(X(:,1), X(:,2)); hold on
plot(L(:,1), L(:,2), 's')

axis equal
plot(-mu,0,'o');
plot(1-mu,0,'o');





