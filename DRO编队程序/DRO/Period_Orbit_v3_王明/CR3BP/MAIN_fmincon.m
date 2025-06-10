
function MAIN_fmincon
%
% 三维周期轨道修正
%
% L4收敛解
% xi = [0.0594550329295266
%     0.86603
%     -8.87166924244361e-24
%     -0.232578397230741
%     0.200899005541197
%     3.75883290145844e-24];
% periOrb_P = 6.55460826598705
% 
% 参考：
% [1]: Vaquero, SpaceCraft Transfer Trajectory Design Exploiting Resonant Orbits in
% Multi-Body Environments [D], Purdue University, 2013, page 33-35
% 
% 作者：张晨
% 2021年5月30日
%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('tools');

%% read para
aux = readInputFile([]);

% 构造初值
if strcmp(aux.orb_IO , 'periOrb')
    [xj_guess , aux] = getIC_periOrb(aux);
elseif strcmp(aux.orb_IO , 'resoOrb')
    [xj_guess , aux] = getIC_resoOrb(aux);
else
    frpintf('wrong!');
end

% 求解NLP
aux.opts = optimoptions(@fmincon , 'Algorithm' , aux.nlp.algorithm , 'Display' , 'iter' , ...
    'TolX' , aux.nlp.tolX , 'TolCon' , aux.nlp.tolCon , 'TolFun' , aux.nlp.tolFun , ...
    'MaxIter' , aux.nlp.maxIter , 'MaxFunEvals' , aux.nlp.maxFunEvals , ...
    'GradObj' , 'on' , 'GradConstr' , 'on' , 'DerivativeCheck' , 'off');
aux.pseudoArc_IO = 0;
[xj , fobj , ExitFlag] = fmincon(@objfcn , xj_guess, [], [], [], [], [], [], @constr, aux.opts, aux);

% -------------------------- 画图 --------------------------
[cc, cc_eq, Gc, Gc_eq] = constr (xj, aux);
plot_solution(xj , Gc_eq , aux);

% -------------------------- 保存结果 --------------------------
save xj xj
save Gc_eq Gc_eq


end
