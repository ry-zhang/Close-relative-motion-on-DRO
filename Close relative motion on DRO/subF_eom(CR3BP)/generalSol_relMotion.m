function rel_motion = generalSol_relMotion(t_sample,k0,k1,k2,theta1_0,theta2_0,para,coe)
%% caculate the general bounded relative motion on DRO
% 
% rel_motion = general_sol_relMotion(t_sample,k0,k1,k2,theta1_0,theta2_0,para,coe)
%
% Input arguments:
% -------------------------------------------------------------------------
% t_sample              [1xn]       要求的相对运动时刻
% k0                    [1x1]       周期分量振幅
% k1                    [1x1]       平面拟周期分量振幅
% k2                    [1x1]       法向拟周期分量振幅
% theta1_0              [1x1]       平面拟周期初始相位
% theta2_0              [1x1]       法向拟周期初始相位
% para.T0/T1/T2         [1x1]       DRO/平面拟周期/法向拟周期 周期
% coe.N                 [1x1]       傅里叶变换采样数
% coe.c3_e1hat          [Nx6]       周期分量傅里叶变换系数矩阵
% coe.c1_e1hat/c1_e2hat [Nx6]       平面拟周期分量傅里叶变换系数矩阵
% coe.c5_e1hat/c6_e2hat [Nx6]       法向拟周期分量傅里叶变换系数矩阵
% 
% Output arguments:
% -------------------------------------------------------------------------
% rel_motion            [6xn]       输出的相对运动状态（归一化单位）
% 
% External functions called:
% -------------------------------------------------------------------------
% iDFTmatrix_theta
% 
% Copyright (C) 1/9/2021 by Chihang Yang 
% email: ychhtl@foxmail.com

%% code

% periodic component
if k0 ~= 0
    e3_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e3hat)';
    rv_rel_periodic = k0*e3_hat_refit;
else
    rv_rel_periodic = zeros(6,length(t_sample));
end

% planar quasi-periodic component
if k1 ~=0
    e1_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e1hat)';
    e2_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e2hat)';
    rv_rel_planar = k1*(cos(theta1_0+t_sample/para.T1*2*pi).*e1_hat_refit + sin(theta1_0+t_sample/para.T1*2*pi).*e2_hat_refit);
else
    rv_rel_planar = zeros(6,length(t_sample));
end

% vertical quasi-periodic component
if k2 ~=0
    e5_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e5hat)';
    e6_hat_refit = real(iDFTmatrix_theta(coe.N,t_sample*2*pi/para.T0) * coe.c1_e6hat)';
    rv_rel_vertical = k2*(cos(theta2_0+t_sample/para.T2*2*pi).*e5_hat_refit + sin(theta2_0+t_sample/para.T2*2*pi).*e6_hat_refit);
else
    rv_rel_vertical = zeros(6,length(t_sample));
end

rel_motion = rv_rel_periodic + rv_rel_planar + rv_rel_vertical;