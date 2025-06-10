clear; clc;
% find periodic orbit

mu = 0.01215; % 20200531

opts_fsolve = optimoptions('fsolve','Display','iter','Algorithm','trust-region',...
    'FunctionTolerance',1e-18,'StepTolerance',1e-13,'MaxIterations',30);
opts1 = odeset('RelTol',1e-13,'AbsTol',1e-20,'Events',@secYeq0);

%% 给出初值
% for ijk = [3 6 11 19 37 98] % 拟周期轨道
load data_DRO_
% 选取几组延拓出的解，给出其周期比、能量、x0
data_sample = 1:10:1201;
T_DRO_sample1 = dataDRO(data_sample,7)/2/pi;
J_sample1 = dataDRO(data_sample,8);
x0_sample1 = dataDRO(data_sample,1);
T_DRO_sample2 = [1/15,1/19,1/24,1/27,1/32,1/36,1/42,1/48,1/54,0.1/27.2844292118811]';
J_sample2 = [3.21651939475916,3.27551907403493,...
    3.34275818065234,3.38050736236399,3.44007802533668,3.48526777267907,...
    3.54973469805399,3.61092835621609,3.66942315218659,5.15310928005317]';
x0_sample2 = [0.948542566656328,0.954526317206957,...
    0.959506515750577,0.961716353812614,0.964593455934316,0.966393866684911,...
    0.968535308630283,0.970212270585268,0.971567330861690,0.982371891725129]';
T_DRO_sample = [T_DRO_sample1;T_DRO_sample2]; % 2~1/54
J_sample = [J_sample1;J_sample2];
x0_sample = [x0_sample1;x0_sample2];
% T_DRO_sample = [5/6,3/4,2/3,1/2,pi/8,1/3,1/4,1/5,1/6,1/7,1/8,
%     1/9,1/10,1/12,1/15,1/19,1/24,1/27,1/32,1/36,1/42,1/48,1/54];
% J_sample = [2.84275066735312,2.87523540955138,2.89674215064624,2.9305237452089,...
%     2.95394755678542,2.97001118561686,3.00048914754382,3.02705366770963,...
%     3.05111699590076,3.07336209842591,3.09419982316304,3.11390471752971,...
%     3.13267202467340,3.16793959925793,3.21651939475916,3.27551907403493,...
%     3.34275818065234,3.38050736236399,3.44007802533668,3.48526777267907,...
%     3.54973469805399,3.61092835621609,3.66942315218659];
% x0_sample = [0.660506774808326,0.712619096702657,0.750480658274125,0.808936204314295,...
%     0.841857960997863,0.85972976187477,0.885197590838893,0.901042109531437,...
%     0.912000610844629,0.920096216952606,0.926356184119068,0.931362405645980,...
%     0.935470965026371,0.941844648311061,0.948542566656328,0.954526317206957,...
%     0.959506515750577,0.961716353812614,0.964593455934316,0.966393866684911,...
%     0.968535308630283,0.970212270585268,0.971567330861690];
% 对x0和能量进行多项式拟合，如此，可以给定周期比，给出可行的x0与能量初值
pJ = fit(T_DRO_sample,J_sample,'cubicinterp');
px = fit(T_DRO_sample,x0_sample,'cubicinterp');
% px = polyfit(T_DRO_sample,x0_sample,6);
% plot(pJ,T_DRO_sample,J_sample)
% 选取目标DRO轨道的周期比（T_DRO/T_Moon）
% 这里为批量优化，选择了一批目标能量，可以根据
% num = 100;
% T_DRO_all = linspace(min(T_DRO_sample),max(T_DRO_sample),num);
% T_DRO_all = linspace(T_DRO_all(ijk),T_DRO_all(ijk+1),num);
% T_DRO_all = linspace(T_DRO_all(3),0.141,num);

T_Moon = 27.2844292118811;
T_MoonSun_synodic = 29.5;
% T_DRO_ratio_str = '[3/4,2/3,3/5,1/2,2/5,1/3,1/4,1/5,1/6]';
% T_DRO_ratio = str2num(T_DRO_ratio_str);

T_DRO_ratio = linspace(27.2844292118811/27.2844292118811,0.1/27.2844292118811,300);
% T_DRO_ratio = linspace(27.2844292118811/27.2844292118811,27/27.2844292118811,300);
% T_DRO_ratio = linspace(27/27.2844292118811,0.1/27.2844292118811,300);
% T_DRO_ratio = [10.60456,7.11394,5.48864]/27.2844292118811;
% T_DRO_ratio = linspace(27.28,0.1,100)/27.2844292118811;
% T_DRO_ratio = linspace(1/54,0.1/27.2844292118811,50);
% T_DRO_ratio = log10(linspace(100^(20/27.2844292118811),100^(1/27.2844292118811),100))/log10(100); % 为了画DRO的图
% T_DRO_all = T_MoonSun_synodic/T_Moon*T_DRO_ratio;
T_DRO_all = T_DRO_ratio;
num = length(T_DRO_all);

J_all = feval(pJ,T_DRO_all);
x0_all = feval(px,T_DRO_all);

% 优化
J_period_all = zeros(num,2);
state_ini_all = zeros(num,4);
delta_T_all = zeros(num,1);

parfor jj = 1:num
    
    disp(num2str(jj))
    
    T_DRO_desired = T_DRO_all(jj);
    J_Initial = J_all(jj);
    x0 = x0_all(jj);

    J_final = fsolve(@(x)continuation(x,x0,T_DRO_desired,opts1,mu),J_Initial,opts_fsolve);

    [delta_T,J_period,state_ini] = continuation(J_final,x0,T_DRO_desired,opts1,mu);
    
    delta_T_all(jj) = delta_T;
    J_period_all(jj,:) = J_period;
    state_ini_all(jj,:) = state_ini;
    
end
% save('DRO_all', 'J_period_all', 'state_ini_all', 'delta_T_all')
% save(['DRO_all_',num2str(ijk)], 'J_period_all', 'state_ini_all', 'delta_T_all')
save('DRO_all44', 'J_period_all', 'state_ini_all', 'delta_T_all',...
    'T_DRO_all')
% save('DRO_all_SEsyPeriod', 'J_period_all', 'state_ini_all', 'delta_T_all',...
%     'T_DRO_all','T_DRO_ratio_str')
% save('DRO_all_all', 'J_period_all', 'state_ini_all', 'delta_T_all',...
%     'T_DRO_all')
% x0_DRO_all = state_ini_all;
% C_DRO_all = J_period_all(:,1);
% T_DRO_all = J_period_all(:,2)*2*pi;
% save('DRO_all3', 'x0_DRO_all', 'C_DRO_all', 'T_DRO_all')





