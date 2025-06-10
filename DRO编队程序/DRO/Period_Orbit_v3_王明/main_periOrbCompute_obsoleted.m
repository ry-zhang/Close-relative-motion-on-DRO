clc;clear;close all
% 三维周期轨道修正
%
% ---------------------- L4周期解 ----------------------
% xi = [0.0594550329295266
%     0.86603
%     -8.87166924244361e-24
%     -0.232578397230741
%     0.200899005541197
%     3.75883290145844e-24];
% periOrb_P = 6.55460826598705
% aux.delta_s = 10e-2;
% nLoop = 100
% ---------------------- 3:2 RSO ----------------------
% 周期轨道初值：
% periOrb_xi =
%     0.3909
%    -0.0000
%     0.0000
%    -0.0000
%     1.5009
%     0.0000
% 周期轨道周期：
% periOrb_P =
%    12.5751
% aux.delta_s = 2e-2;
% nLoop = 1000
% N  = 20
% --------------------
%------------------------2:1DRO-----------------------
% periOrb_xi = 0.808932810866206
%                      0	
%                      0	
%                      0	
%                      0.515640610508855	
%                       0
% periOrb_P =
%    1.57079632679490
% 参考：
% [1]: Vaquero, SpaceCraft Transfer Trajectory Design Exploiting Resonant Orbits in
% Multi-Body Environments [D], Purdue University, 2013, page 33-35
%
% 作者：张晨
% 2021年5月30日
%%%%%%%%%%%%%%%%%%%%%%%%

set(0, 'DefaultAxesFontSize', 14)
set(0, 'DefaultTextFontSize', 14)
set(0, 'DefaultLineLineWidth', 2)

dbstop if error
warning('off');
close all

clc; clear; close all;
addpath('tools');
addpath('CR3BP')
%% read para
aux = readInputFile([]);
if aux.err_IO == 1
    return;
end
%% 设定开关函数初值
aux.pseudoArc_IO = 0;
aux.UseSymmetric_IO = 0;                                    %使用轨道对称性进行校正
aux.Halopseudo = 0;                                               % 专门控制Halo生成
% 构造初值
 [xj , aux] = getIC_periOrb(aux);

fprintf(' ------- 多步打靶 + 最小二乘 --------- \n')

err = inf;
k = 0;
while err > 1e-10  
    aux.pseudoArc_IO = 0;
    
    [~, cc_eq, ~, Gc_eq] = constr(xj, aux);
    Gc_eq = Gc_eq';
    
    % 计算误差
    err = max(abs(cc_eq));
    fprintf('err: %0.4e \n' , err)
    
    % 更新xj
    xj = xj - Gc_eq' * pinv(Gc_eq * Gc_eq') * cc_eq; 
    k = k+1;
end

%% -------------------------- 画图 --------------------------
plot_solution(xj , Gc_eq , aux);

%% --------------------- 伪弧长修正 ----------------------
aux.pseudoArc_IO = 1;
% 赋初值
numContDrect = 1; % 延拓方向
aux.delta_s = 0.1;
if strcmp(aux.orb_type , 'TriangleOrb')
    numContDrect = 1; % 延拓方向
    aux.delta_s = 0.08;
end
if strcmp(aux.orb_type , 'DRO')
    numContDrect = 1; % 延拓方向
    aux.delta_s = 0.1;
end
if strcmp(aux.orb_type , 'Halo')
    numContDrect = -1; % 延拓方向
    if  strcmp(aux.Halo_type , 'L1')
        aux.delta_s = 0.01;
    else
        aux.delta_s = 0.015;
    end
end
if strcmp(aux.orb_type , 'Lyapunov')
    numContDrect = -1; 
    aux.delta_s = 0.05;
     if  strcmp(aux.Lyapunov_type , 'L1')
          
        aux.delta_s = 0.06;
    end
    if  strcmp(aux.Lyapunov_type , 'L3')
     numContDrect = -1; 
        aux.delta_s = 0.08;
    end
end
if strcmp(aux.orb_type , 'ResoOrb')
    numContDrect = 1; % 延拓方向
    aux.delta_s = 0.08;
end
if  strcmp(aux.orb_type , 'LoPO')
    numContDrect = -1; % 延拓方向
    aux.delta_s = 0.05;
end

aux.delta_s  = aux.delta_s *numContDrect;
nLoop = 200;
aux.computeBurf = 0; % 1 = caculating burf
aux.svdBurfTol = 1.4e-3; % tolerance below which impiles a burfication


currBurf = 0; % locally only burf once, donot modify

[~,~,V] = svd(Gc_eq);

aux.xj = xj; % 优化变量
aux.delta_xj = V(:, end); % 零空间向量
% dataSave = zeros(nLoop, 9);
[locs,DRO] = Origin_Orbit_family(120 , aux);
aux.Burf_3D = 0;
% 分岔到三维轨道
if aux.Burf_3D == 1
    if strcmp(aux.orb_type , 'Halo')
        fprintf('--------------------------------------------------- \n');
        fprintf('Halo已经是三维轨道，不必计算分岔 \n');
        return;
    end

    if ~isempty(locs)
        aux.computeBurf = 1; % 1 = caculating burf
        aux.svdBurfTol = 1.4e-3; % tolerance below which impiles a burfication
        
        %% 注意：Lyapunov L1/L2给0.02； L3 给 0.1，DRO,DPO给0.1, LoPO给0.05（Lyapunov），共振给0.1
        aux.burf_delta_s = -0.1; % burf continuation step size
        for ii = 1:length( locs)
            aux.burfLoc = locs(ii);
            dataSave_3D = Burf3D_Orbit( locs(ii)+200,currBurf , aux);
          
        end
%        dataSave_3D2 = Burf3D_Orbit_(locs(ii)+nLoop,currBurf , aux);
    end
end

