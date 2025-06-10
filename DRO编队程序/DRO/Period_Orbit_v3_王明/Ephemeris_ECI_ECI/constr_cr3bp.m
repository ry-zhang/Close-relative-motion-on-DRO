function [cc, cc_eq, Gc, Gc_eq] = constr_cr3bp (yy , aux)
%
% Constraints vector that considers the final point of the trajectory
% (xxn) as a NLP variable (i.e., there is no final integration to impose
% the right boundary condition).
%
% 作者：张晨
% 2021年5月30日
%%%%%%%%%%%%%%%%%%%%%%%%

cc = [];
cc_eq = [];
Gc = [];
Gc_eq = [];

% 拼接点数
n = aux.periOrb_n-1;
if n<2
    printf('n必须大于2')
    return
end
% 变量数
nVar = size(yy , 1);

% 等式约束数

nCeq = 6 * n;


cc_eq = zeros(nCeq , 1);

% 提取周期
periOrb_P = aux.Period;

% reconstruction of the time intervals vector
ttVec = linspace(0 , periOrb_P , n + 1);

% Jacobian of constraints
if nargout > 2
    Gc_eq = zeros(nCeq , nVar);
end

% constraints for the intermediate arcs
for k = 1 : n
    
    % tk & tk1
    tk = ttVec(k);
    tkp1 = ttVec(k + 1);
    
    % x(tk)
    xxk = yy(6 * (k - 1) + 1 : 6 * k);
    
    % varphi(xk , tk : tk1)
    phi0 = eye(6);
    options = odeset('Reltol' , aux.tol , 'AbsTol' , aux.tol);
    [tt_rot , xx_rot] = ode113(@crtbpEqmSTM3D , [tk , tkp1] , [xxk ; phi0(:)] , options , aux);
    
    % 【测试】
    %     plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'b');
    %     plot3(xx_rot(1 , 1) , xx_rot(1 , 2) , xx_rot(1 , 3) , 'bo');
    %     plot3(xx(end , 1) , xx(end , 2) , xx(end , 3) , 'bv');
    
    % STM
    xxkp1_ = xx_rot(end , 1 : 6)';
    phi_tk_tkp1 = reshape(xx_rot(end , 7 : 42) , 6 , 6);
    
    % 第n个节点的约束和雅克比
    if k == n
        
        cc_eq(6 * (k - 1) + 1 : 6 * k) = [xxkp1_(1) - yy(1);
            xxkp1_(2) - yy(2);
            xxkp1_(3) - yy(3);
            xxkp1_(4) - yy(4);
            yy(2) - aux.yFix;
            xxkp1_(6) - yy(6)];
        
        % compute equality constraints Jacobian
        if nargout > 2
            
            % defect constraints Jacobian wrt xxk and xxkp1
            Gc_eq(6*(k-1)+1: 6*k, 6*(k-1)+1: 6*k) = phi_tk_tkp1;
            Gc_eq(6*(k-1)+1: 6*k, 1 : 6) = - eye(6);
            
%             % defect constraints Jacobian wrt P
%             f_xxk = crtbpEqm3D(tk , xxk , aux);
%             f_xxkp1_ = crtbpEqm3D(tkp1 , xxkp1_ , aux);
%             
%             Gc_eq(6*(k-1)+1: 6*k , 6 * n + 1) = -phi_tk_tkp1 * f_xxk * ((k - 1) / n) + f_xxkp1_ * (k / n);
            
            % modi
            Gc_eq(6*(k-1)+5 , :) = 0;
            Gc_eq(6*(k-1)+5 , 2) = 1;
%             Gc_eq(6*(k-1)+5 , 6 * n + 1) = 0;
            
        end
        
    else % 1 : n-1个节点的约束和雅克比
        
        % defect constraints
        % x(tk1)
        xxkp1 = yy(6 * k + 1 : 6 * k + 6);
        cc_eq(6 * (k - 1) + 1 : 6 * k) = xxkp1_ - xxkp1;
        
        % compute equality constraints Jacobian
        if nargout > 2
            
            % defect constraints Jacobian wrt xxk and xxkp1
            Gc_eq(6*(k-1)+1: 6*k, 6*(k-1)+1: 6*k) = phi_tk_tkp1;
            Gc_eq(6*(k-1)+1: 6*k, 6*k+1: 6*(k+1)) = - eye(6);
            
            % defect constraints Jacobian wrt P
%             f_xxk = crtbpEqm3D(tk , xxk , aux);
%             f_xxkp1_ = crtbpEqm3D(tkp1 , xxkp1_ , aux);
%             
%             Gc_eq(6*(k-1)+1: 6*k , 6 * n + 1) = -phi_tk_tkp1 * f_xxk * ((k - 1) / n) + f_xxkp1_ * (k / n);
            
        end
        
    end
    
end

% % 【测试】
% figure(2);
% spy(Gc_eq)

% 等式约束数
% if aux.pseudoArc_IO == 1
%     
%     cc_eq(end , :) = (yy - aux.xj)' * aux.delta_xj - aux.delta_s;
%     Gc_eq(end , :) = aux.delta_xj;
%     
% end

Gc_eq = Gc_eq';
Gc = Gc';

end
