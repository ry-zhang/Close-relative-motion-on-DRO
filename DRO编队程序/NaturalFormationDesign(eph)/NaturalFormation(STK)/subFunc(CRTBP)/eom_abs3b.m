function ydot = eom_abs3b(t,y,mu)
% equation of motion of CR3BP in synodic coordiante system cetered in Moon
%
% ydot = eom_abs3b(t,y,mu)
%
% Input arguments:
% -------------------------------------------------------------------------
% t          [1x1]               FFT size
% y          [6x1]               Position and velocity
% mu         [1x1]               
% 
% Output arguments:
% -------------------------------------------------------------------------
% ydot       [6x1]        dot Matrix
% 
% External functions called:
% -------------------------------------------------------------------------
%  none
% 
% Copyright (C) 25/7/2020 by Chihang Yang 
% email: ychhtl@foxmail.com
% v2(2022-03-01):
% 将旋转坐标系由地球在右边修改为地球在左边

%% 主星/标称轨道在月球会合坐标系M下的动力学
% 原点位于月球旋转坐标系M，x轴从月球指向地球，y轴地球绕月球旋转的方向
r_C = y(1:3); % 主星/标称轨道在月球旋转坐标系M下的位置速度
v_C = y(4:6);
r_C_dot = v_C;
% 圆型限制性三体下
 % 到地球及月球距离取决于会合坐标系的原点是放在质心还是放在月球等
r_E_C = r_C-[-1;0;0]; r_E_C_norm  = norm(r_E_C); r_E_C3 = r_E_C_norm^3;
r_M_C_norm  = norm(r_C); r_M_C3 = r_M_C_norm^3;

Ar_C = [1-mu/r_M_C3-(1-mu)/r_E_C3, 0, 0;
        0, 1-mu/r_M_C3-(1-mu)/r_E_C3, 0;
        0, 0, -mu/r_M_C3-(1-mu)/r_E_C3];
Av_C = [0, 2, 0;
        -2, 0, 0;
        0, 0, 0];
B = [(1-mu)*(1-1/r_E_C3); 0; 0];
v_C_dot = Ar_C*r_C + Av_C*v_C + B;
ydot = [r_C_dot; v_C_dot];
% 
% M1 = -1/r_M_C_norm^5 * (-3*r_C*r_C'+eye(3)*r_M_C_norm^2);
% M2 = -1/r_E_C_norm^5 * (-3*r_E_C*r_E_C'+eye(3)*r_E_C_norm^2);
% v_C_ddot = Ar_C*v_C + (Av_C - mu*M1 - (1-mu)*M2)*v_C_dot;