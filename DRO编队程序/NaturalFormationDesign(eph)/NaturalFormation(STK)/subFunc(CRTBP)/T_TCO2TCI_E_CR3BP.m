function [rv_TCI_c] = T_TCO2TCI_E_CR3BP(rv_TCO_c,rv_ECI_t,t_sample,Flag_frame,mu)
% 主星原点惯性系(target-centered inertial,TCO)转化到轨道坐标系(target-centered orbital,TCI)
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH
% v3 2021/8/20 YCH

%输入：
%    跟踪航天器相对位置在ECI中的分量列阵rv_ECI_t
%    目标航天器位置速度rMT_ECI,vMT_ECI
%    坐标系标志，('VNC','LVLH','VVLH')
%输出：
%    跟踪航天器相对位置在TC中的分量列阵rv_TCO_c

%轨道坐标系VNC
%    -原点位于目标航天器质心
%    -x轴沿着目标航天器的速度方向
%    -y轴沿着目标航天器相对于中心天体的角速度方向 
%轨道坐标系LVLH(即RIC)
%    -原点位于目标航天器质心
%    -x轴沿着中心天体指向目标航天器
%    -z轴沿着目标航天器相对于中心天体的角速度方向
%轨道坐标系VVLH
%    -原点位于目标航天器质心
%    -y轴沿着目标航天器相对于中心天体的角速度的负方向
%    -z轴沿着目标航天器指向中心天体

%地心惯性系ECI
%    -原点位于地心
%    -x轴为0时刻的地月连线方向（地指向月）
%    -z轴为月球相对于地球的角动量方向

% 判断输入的rv_TCO_c是否包含速度
[row,col] = size(rv_TCO_c);
if col == 3
    rv_TCI_c = zeros(row,3);
elseif col == 6
    rv_TCI_c = zeros(row,6);
else
    error('Wrong input size')
end

for ii = 1:row
    r_TCO_c_ii = rv_TCO_c(ii,1:3); 
    if col == 3
        v_TCO_c_ii = [];
    else
        v_TCO_c_ii = rv_TCO_c(ii,4:6);
    end
    r_ECI_t_ii = rv_ECI_t(ii,1:3); v_ECI_ii = rv_ECI_t(ii,4:6);
    theta = t_sample(ii);
    [r_TC_c_ii,v_TC_c_ii] = TCO2TCI(r_TCO_c_ii',v_TCO_c_ii',r_ECI_t_ii',v_ECI_ii',theta,Flag_frame,mu);
    rv_TCI_c(ii,:) = [r_TC_c_ii',v_TC_c_ii'];
end

end

function [r_TCI_c,v_TCI_c] = TCO2TCI(r_TCO_c,v_TCO_c,r_ECI_t,v_ECI_t,theta,Flag_frame,mu)

% 动力学求vdot
% 圆型限制性三体下
% 地心坐标系，原点位于地心
r_E_t_norm  = norm(r_ECI_t); r_E_t3 = r_E_t_norm^3;
r_M_t = r_ECI_t+[cos(theta);sin(theta);0];
r_M_t_norm  = norm(r_M_t); r_M_t3 = r_M_t_norm^3;

v_ECI_t_dot = -(1-mu)/r_E_t3*r_ECI_t-mu/r_M_t3*r_M_t;

v_ECI_t_norm = norm(v_ECI_t);
h_ECI_t = cross(r_ECI_t,v_ECI_t); 
h_ECI_t_norm = norm(h_ECI_t);
% h_ECI_t = [0,0,1]';
% h_ECI_t_norm = 1;


% 轨道坐标系的单位向量及其导数在旋转坐标系中的表示
if strcmp(Flag_frame,'VNC')
    i_OinR = v_ECI_t/v_ECI_t_norm;
    j_OinR = h_ECI_t/h_ECI_t_norm;
    k_OinR = cross(i_OinR,j_OinR);
    i_OinR_dot = - dot(v_ECI_t,v_ECI_t_dot) / v_ECI_t_norm^3 * v_ECI_t + v_ECI_t_dot/v_ECI_t_norm;
    j_OinR_dot = - dot(v_ECI_t_dot,j_OinR) * r_E_t_norm / h_ECI_t_norm * cross(j_OinR,r_ECI_t/r_E_t_norm);
%     j_OinR_dot = [0,0,0];
    k_OinR_dot = cross(i_OinR_dot,j_OinR)+cross(i_OinR,j_OinR_dot);
%     omegaT = (dot(v_ECI_t_dot,h_ECI_t)*cross(v_ECI_t,h_ECI_t)...
%         + dot(v_ECI_t_dot,h_ECI_t)*dot(r_ECI_t,v_ECI_t)*v_ECI_t...
%         - dot(v_ECI_t_dot,cross(v_ECI_t,h_ECI_t))*h_ECI_t)/h_ECI_t_norm^2/v_ECI_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    k_OinR = h_ECI_t/h_ECI_t_norm;
    i_OinR = r_ECI_t/r_E_t_norm;
    j_OinR = cross(k_OinR,i_OinR);
    i_OinR_dot = dot(v_ECI_t,j_OinR) / r_E_t_norm * j_OinR;
    k_OinR_dot = - dot(v_ECI_t_dot,k_OinR) * r_E_t_norm / h_ECI_t_norm * j_OinR;
%     k_OinR_dot = [0,0,0];
    j_OinR_dot = r_E_t_norm/h_ECI_t_norm * dot(v_ECI_t_dot,k_OinR)*k_OinR - 1/r_E_t_norm * dot(v_ECI_t,j_OinR)*i_OinR;
    % j_OinR_dot2 = cross(k_OinR_dot,i_OinR)+cross(k_OinR,i_OinR_dot);
%     omegaT = h_ECI_t/r_ECI_t_norm^2 + dot(v_ECI_t_dot,h_ECI_t)/h_ECI_t_norm^2*r_ECI_t;
elseif strcmp(Flag_frame,'VVLH')
    j_OinR = -h_ECI_t/h_ECI_t_norm;
    k_OinR = -r_ECI_t/r_E_t_norm;
    i_OinR = cross(j_OinR,k_OinR);
    k_OinR_dot = - dot(v_ECI_t,i_OinR) / r_E_t_norm * i_OinR;
    j_OinR_dot = - dot(v_ECI_t_dot,j_OinR) * r_E_t_norm / h_ECI_t_norm * i_OinR;
%     j_OinR_dot = [0,0,0];
    i_OinR_dot = r_E_t_norm/h_ECI_t_norm * dot(v_ECI_t_dot,j_OinR)*j_OinR + 1/r_E_t_norm * dot(v_ECI_t,i_OinR)*k_OinR;
%     omegaT = h_ECI_t/r_ECI_t_norm^2 + dot(v_ECI_t_dot,h_ECI_t)/h_ECI_t_norm^2*r_ECI_t;
else
    error('Wrong orbital frame flag')
end

%坐标转换
TM_TCO2TCI = [i_OinR j_OinR k_OinR];
r_TCI_c = TM_TCO2TCI*r_TCO_c;

if isempty(v_TCO_c)
    v_TCI_c = [];
else
    TM_TCO2TCI_dot = [i_OinR_dot,j_OinR_dot,k_OinR_dot];
    v_TCI_c = TM_TCO2TCI*v_TCO_c + TM_TCO2TCI_dot*r_TCO_c;
%     v_TCI_c = TM_TCO2TCI*v_TCO_c + cross(omegaT,r_TCI_c);
end

end

