function [x0,F_final,Iter_error] = fsolve_GN(funname,x_ini,paraIter,varargin)

% Gauss-Newton迭代算法
% [x0,Iter_error,F_final] = fsolve_GN(@funname,x_ini,paraIter,para,paraCT);
%
% Input arguments:
% -------------------------------------------------------------------------
% funname    [function handle]              优化函数名
%            @funname(x0,varargin{:})
% x_ini      [Nx1]                          优化初值
% 
% paraIter   strcture. MaxError / MinDx /   GN迭代设置
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
x0 = x_ini; ii_Newton = 1;
Iter_error = 0.5; dx = paraIter.MinDx+1; 

if ~isfield(paraIter,'Display')
    paraIter.Display = 1;
end

if paraIter.Display == 1
    disp(['  Iteration' , '    MaxError', '     ComputationTime ', '    Cond_H']);
end

%% 迭代
%         while(ii_Newton<MaxIter && Iter_error2>MaxError )
while (ii_Newton<paraIter.MaxIter) && (Iter_error>paraIter.MaxError) && (norm(dx,inf)>paraIter.MinDx)
    tic
    [F,J] = funname(x0,varargin{:});

    % 删除梯度矩阵里全为0的行, 精简梯度矩阵，减小计算量及矩阵的条件数
    JF_filter_row = sum(J,2)~=0;
    J1 = J(JF_filter_row,:);
    F1 = F(JF_filter_row,:);

    JF_filter_col = sum(J1,1)~=0;
    J1 = J1(:,JF_filter_col);
    
    % 牛顿法，分别用GPU和CPU求逆
%     if N>700
%         dx = pinv_GPU(J1)*F;
%     else
%         dx = pinv(J1)*F1;
%         cond_H = 1; 
%     end

    % 拟牛顿法1 - 不缩放矩阵
%     H = J1'*J1;
%     dx = H\(J1'*F1);

    % 拟牛顿法2 - 缩放矩阵以改善条件数
    [row_J1,col_J1] = size(J1);
    key_R = 1./mean(abs(J1),1);  DR = sparse(1:col_J1,1:col_J1,key_R);
    J2 = J1*DR;
    key_L = 1./mean(abs(J2),2);  DL = sparse(1:row_J1,1:row_J1,key_L);
    J3 = DL*J2;
    F3 = DL*F1;
    H = J3'*J3;
    cond_H = condest(H); % H矩阵的条件数，越大则求逆精度越低
    % 当条件数过大时，采用伪逆，条件数较小时，采用高斯牛顿
    if cond_H>1e15
        dx = pinv(J1)*F1;
    else
        dx = key_R'.*(H\(J3'*F3));
    end
    
    try
        x0(JF_filter_col) = x0(JF_filter_col) - dx;
    catch
        x0(JF_filter_col) = x0(JF_filter_col) - dx';
    end

    Iter_error = max(abs(F));
    tCompu = toc;
    if paraIter.Display == 1
        disp(['      ',num2str(ii_Newton), '       ',num2str(Iter_error,'%8.4e'),...
            '        ', num2str(tCompu,'%8.4f'),'s','       ',num2str(cond_H,'%8.4e') ]);
    end
    ii_Newton = ii_Newton+1;
    if Iter_error>paraIter.BreakError
        break
    end
end
F_final = F;
