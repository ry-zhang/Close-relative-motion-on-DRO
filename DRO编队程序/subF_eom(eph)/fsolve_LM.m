function [x0,F_final,Iter_error] = fsolve_LM(funname,x_ini,paraIter,varargin)

% Levenberg-Marquardt迭代算法，结合高斯牛顿与梯度下降，
% 收敛性强，在初值离较远时计算速度快
% https://www.cnblogs.com/monoSLAM/p/5249339.html
% [x0,Iter_error,F_final] = fsolve_LM(@funname,x_ini,paraIter,para,paraCT);
%
% Input arguments:
% -------------------------------------------------------------------------
% funname    [function handle]              优化函数名
%            @funname(x0,varargin{:})
% x_ini      [Nx1]                          优化初值
% 
% paraIter   strcture. MaxError / MinDx /   LM迭代设置
%            BreakError / MaxIter
% varargin                                  funname除x0外的其他输入
% 
% Output arguments:
% -------------------------------------------------------------------------
% x0          [Nx1]                         迭代结果
% F_final     [Mx1]                         方程组迭代误差
% Iter_error  [1x1]                         norm(F_final,inf);
% 
% External functions called:
% -------------------------------------------------------------------------
%  funname
% 
% Copyright (C) 2021-03-02 by Chihang Yang 
% email: ychhtl@foxmail.com

%% 迭代初始化
Iter_error = 1; ii_LM = 1;
x0 = x_ini; dx = paraIter.MinDx+1; 
lambda = 0.01;  nu = 2;
updateJ = 1;

if ~isfield(paraIter,'Display')
    paraIter.Display = 1;
end

if paraIter.Display == 1
    disp(['  Iteration', '    MaxError', '       lambda','           g', ...
            '           dx','        ComputationTime', '     Cond_H_LM']); 
end

%% 迭代
while(ii_LM<paraIter.MaxIter) && (Iter_error>paraIter.MaxError) && (norm(dx,inf)>paraIter.MinDx)
    tic
    if updateJ==1
        
        [F,J] = funname(x0,varargin{:});
        
        % 精简梯度矩阵，减小计算量及矩阵的条件数
        % 删除梯度矩阵里全为0的行
        JF_filter_row = sum(J,2)~=0;
        J1 = J(JF_filter_row,:);
        F1 = F(JF_filter_row,:);
        % 删除梯度矩阵里全为0的列
        JF_filter_col = sum(J1,1)~=0;
        J1 = J1(:,JF_filter_col);

        H = J1'*J1; % 拟海塞矩阵
        Iter_error = norm(F,inf);
    end
    H_LM = H+lambda*eye(size(H));
    g = (J1'*F1);
    cond_H_LM = condest(H_LM); % H矩阵的条件数，越大则求逆精度越低
    dx = H_LM\g;

    % 更新x_LM
    x_LM = x0;
    try
        x_LM(JF_filter_col) = x0(JF_filter_col) - dx;
    catch
        x_LM(JF_filter_col) = x0(JF_filter_col) - dx';
    end

    F_LM = funname(x_LM,varargin{:});
    Iter_error_LM = norm(F_LM,inf);

    rho = (norm(F,inf)-norm(F_LM,inf))/(dx'*(lambda*dx+g));
    if Iter_error_LM<Iter_error
        lambda = lambda*max(1/10,1-(2*rho-1)^3);
        nu = 2;
        x0 = x_LM;
        Iter_error = Iter_error_LM;
        updateJ = 1;
        F_final = F_LM;
    else
        lambda = lambda*nu; nu = nu*2;
        updateJ = 0;
        F_final = F;
    end
    tCompu = toc;
    if paraIter.Display == 1
        disp(['      ', num2str(ii_LM),'       ', num2str(Iter_error,'%8.4e'),'    ',...
            num2str(lambda,'%8.4e'), '    ', num2str(norm(g,inf),'%8.4e'), '   ', ...
            num2str(norm(dx),'%8.4e'),'        ', num2str(tCompu,'%8.4f'),'s', ...
            '        ', num2str(cond_H_LM,'%8.4e')]);
    end
    
    ii_LM = ii_LM+1;
    if Iter_error>paraIter.BreakError
        break
    end
end