function [cc, cc_eq, Gc, Gc_eq] = constr (xx, aux)
%
% 周期轨道四体修正
%
% 作者：王明
% 2021年11月4日
%%%%%%%%%%%%%%%%%%%%%%%%

cc = [];
cc_eq = [];
Gc = [];
Gc_eq = [];
% number of nodes

N = aux.node_n;

% cal state & time of flight & epoch
XX_mtx = reshape(xx(1 : 6 * N) , 6 , N);
T_vec = xx(6 * N + 1 : 6 * N + (N - 1));
tau_vec = xx(6 * N + (N - 1) + 1 : 6 * N + (N - 1) + N);

% number of vars
nVar = size(xx , 1);

% number of constraints
nCeq = 7 * (N - 1);
cc_eq = zeros(nCeq , 1);

% Jacobian of constraints
Gc_eq = zeros(nCeq , nVar);

% 保存结果
if aux.plot_IO == 1
    result_tt = [];
    result_yy = [];
end
% constraints for multi-shooting
for iLoop = 1 : (N - 1)
    
    % load vars
    xxk = XX_mtx(: , iLoop);
    xxkp1 = XX_mtx(: , iLoop + 1);
    tauk = tau_vec(iLoop);
    taukp1 = tau_vec(iLoop + 1);
    Tk = T_vec(iLoop);
    phi0 = eye(6);
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
     [tt , yy] = ode113(@eqmJ2000STM , [tauk , tauk + Tk] , [xxk ; phi0(:) ; zeros(6 , 1)] , options , aux);
     if aux.plot_IO == 1
         result_tt = [result_tt ; tt];
         result_yy = [result_yy ; yy];
     end
    % STM
    xxkp1_ = yy(end , 1 : 6)';
    phi_tk_tkp1 = reshape(yy(end , 7 : 42) , 6 , 6);
    
    % dx_dtauk
    dx_dtauk = yy(end , 43 : 48);
    % ---------------- defect constraints ----------------
    % link constraint
    cc_eq(6 * (iLoop - 1) + 1 : 6 * iLoop) = xxkp1_ - xxkp1;
    
    % epoch constraint, tauk, Tk, tauk1
    cc_eq(6 * (N - 1) + iLoop , :) = tauk + Tk - taukp1;
       % compute equality constraints Jacobian
    if nargout > 2
        
        % link constraint wrt xxk
        Gc_eq(6 * (iLoop - 1) + 1 : 6 * iLoop , 6 * (iLoop - 1) + 1 : 6 * iLoop) = phi_tk_tkp1;
        
        % link constraint wrt xxkp1
        Gc_eq(6 * (iLoop - 1) + 1 : 6 * iLoop , 6 * iLoop + 1 : 6 * (iLoop + 1)) = - eye(6);
        
        % link constraint wrt Tk
         f_xxkp1_ = eqmJ2000STM(tauk + Tk , xxkp1_ , aux);
        Gc_eq(6 * (iLoop - 1) + 1 : 6 * iLoop , 6 * N + iLoop) = f_xxkp1_;
        
        % link constraint wrt tauk
        Gc_eq(6 * (iLoop - 1) + 1 : 6 * iLoop , 6 * N + (N - 1) + iLoop) = dx_dtauk;
        
        % epoch constraint wrt Tk
        Gc_eq(6 * (N - 1) + iLoop , 6 * N + iLoop) = 1;
        
        % epoch constraint wrt tauk
        Gc_eq(6 * (N - 1) + iLoop , 6 * N + (N - 1) + iLoop) = 1;
        
        % epoch constraint wrt tauk1
        Gc_eq(6 * (N - 1) + iLoop , 6 * N + (N - 1) + iLoop + 1) = -1;
        
    end
end
% 额外约束需求
if strcmp(aux.method , 'newton')
%         Gc_eq(:,2) = 0;
%         Gc_eq(:,6 * N + (N - 1) + 1) = 0;
%         Gc_eq(:,6*(N-1)+2) = 0;
    x0 = [xx(1:3)*aux.LU;  xx(4:6)*aux.VU];
    xn =  [xx(6*(N-1)+1:6*(N-1)+3)*aux.LU;  xx(6*(N-1)+4:6*(N-1)+6)*aux.VU];
    x0_rot = Eme2BRC(aux.dep_t_jdtdb +tau_vec(1) * aux.TU / 86400, x0, aux);
    T_Mtx = get_Eme2RotJcoMtx(tau_vec(1), aux);
     T_Mtx =  T_Mtx^-1;
     xn_rot = Eme2BRC(aux.dep_t_jdtdb +tau_vec(N) * aux.TU / 86400, xn, aux);
%      cc_eq2 = [x0_rot(2)-aux.yFix; xx(6 * N + (N - 1) + 1)];   % 使得y为0，tau1为0
     cc_eq2 = [x0_rot(2)-aux.yFix; xx(6 * N + (N - 1) + 1);xn_rot(2)-aux.yFix];                   % 使得y为0，tau1为0，末端节点y为0.
     Gc_eq2 = zeros(size(cc_eq2,1), nVar);
     Gc_eq2(1,1:6) = T_Mtx(2,:);
     Gc_eq2(2,6 * N + (N - 1) + 1) = 1;
     
     T_Mtx = get_Eme2RotJcoMtx(tau_vec(N), aux);
     T_Mtx =  T_Mtx^-1;
     Gc_eq2(3,6*(N-1)+1:6*(N-1)+6) = T_Mtx(2,:);
     
    cc_eq = [cc_eq; cc_eq2];
    Gc_eq = [Gc_eq; Gc_eq2];
end
if strcmp(aux.method , 'fmincon')
    x0 = [xx(1:3)*aux.LU;  xx(4:6)*aux.VU];
   
    x0_rot = Eme2BRC(aux.dep_t_jdtdb +tau_vec(1) * aux.TU / 86400, x0, aux);
    T_Mtx = get_Eme2RotJcoMtx(tau_vec(1), aux);
    T_Mtx =  T_Mtx^-1;
    
    %     cc_eq2 = [x0_rot(2)-aux.yFix; xx(6 * N + (N - 1) + 1)];
    cc_eq2 = [x0_rot(2)-aux.yFix; xx(6 * N + (N - 1) + 1)];                   % 使得y为0，tau1为0，末端节点y为0.
    Gc_eq2 = zeros(2, nVar);
    Gc_eq2(1,1:6) = T_Mtx(2,:);
    Gc_eq2(2,6 * N + (N - 1) + 1) = 1;
    
    cc_eq = [cc_eq; cc_eq2];
    Gc_eq = [Gc_eq; Gc_eq2];
end
Gc_eq = Gc_eq';
Gc = Gc';
% 收敛后画图
% --------------------- 画图 ---------------------
if aux.plot_IO == 1
    
    % 旋转系画图
    yy_rot = zeros(length(result_tt) , 6);
    for jLoop = 1 : length(result_tt)
        % 【G】
          yy_rot(jLoop , :) = Eme2BRC(aux.dep_t_jdtdb + result_tt(jLoop) * aux.TU / 86400,  [result_yy(jLoop , 1:3) * aux.LU, result_yy(jLoop , 4:6) * aux.VU]', aux);
%         yy_rot(jLoop , :) = crtbpEmeToRot(aux.dep_t_jdtdb + result_tt(jLoop) * aux.TU / 86400 , 3 , 10 , 'P1' , ...
%             [result_yy(jLoop , 1:3) * aux.LU, result_yy(jLoop , 4:6) * aux.VU]' , aux);
    end
    
    figure(1);
    plot3(yy_rot(: , 1) , yy_rot(: , 2) , yy_rot(: , 3) , 'r' , 'linewidth' , 2);
    
    % 惯性系画图
    figure(2);
    plot3(result_yy(: , 1) , result_yy(: , 2) , result_yy(: , 3) , 'r' , 'linewidth' , 2);
    
    % 梯度矩阵
%     h3 = figure(3); hold on; grid on;
%     set(h3 , 'position' , [300 , 300 , 500 , 500])
%     spy(Gc_eq');
    
end


end

