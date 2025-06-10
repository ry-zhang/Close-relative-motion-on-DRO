clc;clear;close all
% 三维周期轨道星历修正
%
% 作者：王明
%
% 2021年11月4日
%
%%%%%%%%%%%%%%%%%%%%%%%%
addpath('tools');

% addpath('Orbit_data');
addpath('Ephemeris_ECI_ECI');
set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultTextFontSize', 14)
set(0, 'DefaultLineLineWidth', 2)

dbstop if error
warning('off');
%% read para
aux = readInputFile([]);
if aux.err_IO == 1
    return;
end
aux.method = 'newton';
% % 构造初值
[xj,aux] = getIC_periOrb(aux);



% %------------------四体修正------------------------------
% step1: 牛顿法粗修正
fprintf(' ------- 多步打靶 + 最小二乘 --------- \n')
err = inf;
tic
iterNumMax = 50;
iterNum = 1;
while err > 5e-10 && iterNum < iterNumMax
    
    % 约束 / 约束梯度
    aux.plot_IO = 0;
    [~, cc_eq, ~, Gc_eq] = constr (xj , aux);
    Gc_eq = Gc_eq';
    % 计算误差
    err = max(abs(cc_eq));
    fprintf('err: %0.4e \n' , err)
    % 更新xj
    xj = xj - Gc_eq' * pinv(Gc_eq * Gc_eq') * cc_eq;
    % xj = xj - Gc_eq \ cc_eq;
    
    iterNum = iterNum + 1;
end
toc
aux.plot_IO = 1;
[~, cc_eq, ~, Gc_eq] = constr (xj , aux);
% x0与xn是被当下的地月距离与月球公转速度归一化的数据
x0 = Eme2BRC(aux.dep_t_jdtdb +xj(7*aux.node_n) * aux.TU / 86400, [xj(1:3)*aux.LU;xj(4:6)*aux.VU], aux);
xn = Eme2BRC(aux.dep_t_jdtdb +xj(8*aux.node_n-1) * aux.TU / 86400, [xj(6*aux.node_n-5:6*aux.node_n-3)*aux.LU;xj(6*aux.node_n-2:6*aux.node_n)*aux.VU], aux);
fprintf('收敛轨道初值：%0.15f, %0.15f, %0.15f, %0.15f, %0.15f, %0.15f \n' , x0)
fprintf('收敛轨道始末差：%0.15f, %0.15f, %0.15f, %0.15f, %0.15f, %0.15f \n' , xn-x0)
fprintf('收敛轨道初始时刻：%0.15f \n' , xj(7*aux.node_n))
fprintf('收敛轨道周期：%0.15f.\n' ,xj(8*aux.node_n-1)- xj(7*aux.node_n))
fprintf('--------------------计算完毕---------------------\n')

jd0 = aux.dep_t_jdtdb;
x0j2k = [xj(1:3)*aux.LU;xj(4:6)*aux.VU];
t_total = (xj(8*aux.node_n-1)- xj(7*aux.node_n))*aux.TU;
save DROephMultiRev2024_T33 jd0 x0j2k t_total

% clf(figure(1));

% % step2: fmioncon精修正
%   'TolX' , 1e-9 , 'TolCon' , 1e-9 , 'TolFun' , 1e-3 ,...
% fprintf(' -------fmincon 局部寻优-------- \n')
% aux.method = 'fmincon';
% aux.opts = optimoptions(@fmincon , 'Algorithm' , aux.nlp.algorithm , 'Display' , 'iter' , ...
%     'MaxIter' , aux.nlp.maxIter , 'MaxFunEvals' , aux.nlp.maxFunEvals , ...
%      'TolX' , 1e-10 , 'TolCon' , 1e-10 , 'TolFun' , 1e-4 ,...
%      'GradObj' , 'off' , 'GradConstr' , 'off ' , 'DerivativeCheck' , 'off' , 'UseParallel' , true);
% aux.plot_IO = 0;
% tic
% [xj , fobj , ExitFlag] = fmincon(@objfun , xj , [] , [] , [] , [] , [] , [] , @constr , aux.opts, aux);
% toc
% aux.plot_IO = 1;
% [~, cc_eq, ~, Gc_eq] = constr (xj , aux);
% x0 = Eme2BRC(aux.dep_t_jdtdb +xj(7*aux.node_n) * aux.TU / 86400, [xj(1:3)*aux.LU;xj(4:6)*aux.VU], aux);
% xn = Eme2BRC(aux.dep_t_jdtdb +xj(8*aux.node_n-1) * aux.TU / 86400, [xj(6*aux.node_n-5:6*aux.node_n-3)*aux.LU;xj(6*aux.node_n-2:6*aux.node_n)*aux.VU], aux);
% fprintf('收敛轨道初值：%0.15f, %0.15f, %0.15f, %0.15f, %0.15f, %0.15f \n' , x0)
% fprintf('收敛轨道始末差：%0.15f, %0.15f, %0.15f, %0.15f, %0.15f, %0.15f \n' , xn-x0)
% fprintf('收敛轨道初始时刻：%0.15f \n' , xj(7*aux.node_n))
% fprintf('收敛轨道周期：%0.15f.\n' ,xj(8*aux.node_n-1)- xj(7*aux.node_n))
% % if strcmp(aux.method , 'newton')
%     fprintf(' ------- 多步打靶 + 最小二乘 --------- \n')
%     err = inf;
%     tic
%     iterNumMax = 50;
%     iterNum = 1;
%     while err > 1e-10 && iterNum < iterNumMax
%         
%         % 约束 / 约束梯度
%         aux.plot_IO = 0;
%         [~, cc_eq, ~, Gc_eq] = constr (xj , aux);
%         Gc_eq = Gc_eq';
%         % 计算误差
%         err = max(abs(cc_eq));
%         fprintf('err: %0.4e \n' , err)
%         % 更新xj
%         xj = xj - Gc_eq' * pinv(Gc_eq * Gc_eq') * cc_eq;
%         % xj = xj - Gc_eq \ cc_eq;
%         
%         iterNum = iterNum + 1;
%     end
%     toc
%     
%     aux.plot_IO = 1;
%     [~, cc_eq, ~, Gc_eq] = constr (xj , aux);
% fprintf('收敛轨道初值：%0.15f, %0.15f, %0.15f, %0.15f, %0.15f, %0.15f \n' , xj(1:6))
% fprintf('收敛轨道初始时刻：%0.15f \n' , xj(7*aux.node_n))
% fprintf('收敛轨道周期：%0.15f.\n' ,xj(8*aux.node_n-1)- xj(7*aux.node_n))
%  fprintf('--------------------计算完毕---------------------\n')
% elseif strcmp(aux.method , 'fmincon')
%     % 【测试】检查梯度！
%     aux.opts = optimoptions(@fmincon , 'Algorithm' , aux.nlp.algorithm , 'Display' , 'iter' , ...
%         'MaxIter' , aux.nlp.maxIter , 'MaxFunEvals' , aux.nlp.maxFunEvals , ...
%         'GradObj' , 'off' , 'GradConstr' , 'off' , 'DerivativeCheck' , 'off' , 'UseParallel' , true);
%     aux.plot_IO = 0;
%     tic
%     [xj , fobj , ExitFlag] = fmincon(@objfun , xj , [] , [] , [] , [] , [] , [] , @constr , aux.opts, aux);
%     toc
%     aux.plot_IO = 1;
%     [~, cc_eq, ~, Gc_eq] = constr (xj , aux);
%     
% end
%    
