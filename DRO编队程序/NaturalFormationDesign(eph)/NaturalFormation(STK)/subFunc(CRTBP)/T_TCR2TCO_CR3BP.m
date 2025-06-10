function [rv_TCO_c] = T_TCR2TCO_CR3BP(rv_TCR_c,rv_MCR_t,Flag_frame,mu)
% 主星原点旋转系(target-centered rotational,TCR)转化到轨道坐标系(target-centered orbital,TCO)
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH

%输入：
%    跟踪航天器相对位置在MCR中的分量列阵rv_MCR_t
%    目标航天器位置速度rMT_MCR,vMT_MCR
%    坐标系标志，('VNC','LVLH','VVLH')
%输出：
%    跟踪航天器相对位置在TC中的分量列阵rv_TCR_c

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

%月心旋转系MCR
%    -原点位于月心
%    -x轴沿月球指向地球方向
%    -z轴为瞬时月球相对于地球的角动量方向

% 判断输入的rv_TCR_c是否包含速度
[row,col] = size(rv_TCR_c);
if col == 3
    rv_TCO_c = zeros(row,3);
elseif col == 6
    rv_TCO_c = zeros(row,6);
else
    error('Wrong input size')
end

for ii = 1:row
    r_TCR_c_ii = rv_TCR_c(ii,1:3); 
    if col == 3
        v_TCR_c_ii = [];
    else
        v_TCR_c_ii = rv_TCR_c(ii,4:6);
    end
    r_MCR_t_ii = rv_MCR_t(ii,1:3); v_MCR_t_ii = rv_MCR_t(ii,4:6);
    [r_TCO_c_ii,v_TCO_c_ii] = TCR2TCO(r_TCR_c_ii',v_TCR_c_ii',r_MCR_t_ii',v_MCR_t_ii',Flag_frame,mu);
    rv_TCO_c(ii,:) = [r_TCO_c_ii',v_TCO_c_ii'];
end

end

function [r_TCO_c,v_TCO_c] = TCR2TCO(r_TCR_c,v_TCR_c,r_MCR_t,v_MCR_t,Flag_frame,mu)

% 动力学求vdot
% 圆型限制性三体下
% 到地球及月球距离取决于会合坐标系的原点是放在质心还是放在月球等
r_E_t = r_MCR_t-[-1;0;0]; r_E_t_norm  = norm(r_E_t); r_E_t3 = r_E_t_norm^3;
r_MCR_t_norm  = norm(r_MCR_t); r_MCR_t3 = r_MCR_t_norm^3;

Ar_MCR = [1-mu/r_MCR_t3-(1-mu)/r_E_t3, 0, 0;
        0, 1-mu/r_MCR_t3-(1-mu)/r_E_t3, 0;
        0, 0, -mu/r_MCR_t3-(1-mu)/r_E_t3];
Av_MCR = [0, 2, 0;
        -2, 0, 0;
        0, 0, 0];
B_MCR = [(1-mu)*(1-1/r_E_t3); 0; 0];
v_MCR_t_dot = Ar_MCR*r_MCR_t + Av_MCR*v_MCR_t + B_MCR;

h_MCR_t = cross(r_MCR_t,v_MCR_t); 
h_MCR_t_norm = norm(h_MCR_t);
v_MCR_t_norm = norm(v_MCR_t);

% 轨道坐标系的单位向量及其导数在旋转坐标系中的表示
if strcmp(Flag_frame,'VNC')
    i_OinR = v_MCR_t/v_MCR_t_norm;
    j_OinR = h_MCR_t/h_MCR_t_norm;
    k_OinR = cross(i_OinR,j_OinR);
    i_OinR_dot = - dot(v_MCR_t,v_MCR_t_dot) / v_MCR_t_norm^3 * v_MCR_t + v_MCR_t_dot/v_MCR_t_norm;
    j_OinR_dot = - dot(v_MCR_t_dot,j_OinR) * r_MCR_t_norm / h_MCR_t_norm * cross(j_OinR,r_MCR_t/r_MCR_t_norm);
    k_OinR_dot = cross(i_OinR_dot,j_OinR)+cross(i_OinR,j_OinR_dot);
%     omegaT = (dot(v_MCR_t_dot,h_MCR_t)*cross(v_MCR_t,h_MCR_t)...
%         + dot(v_MCR_t_dot,h_MCR_t)*dot(r_MCR_t,v_MCR_t)*v_MCR_t...
%         - dot(v_MCR_t_dot,cross(v_MCR_t,h_MCR_t))*h_MCR_t)/h_MCR_t_norm^2/v_MCR_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    k_OinR = h_MCR_t/h_MCR_t_norm;
    i_OinR = r_MCR_t/r_MCR_t_norm;
    j_OinR = cross(k_OinR,i_OinR);
    i_OinR_dot = dot(v_MCR_t,j_OinR) / r_MCR_t_norm * j_OinR;
    k_OinR_dot = - dot(v_MCR_t_dot,k_OinR) * r_MCR_t_norm / h_MCR_t_norm * j_OinR;
    j_OinR_dot = r_MCR_t_norm/h_MCR_t_norm * dot(v_MCR_t_dot,k_OinR)*k_OinR - 1/r_MCR_t_norm * dot(v_MCR_t,j_OinR)*i_OinR;
    % j_OinR_dot2 = cross(k_OinR_dot,i_OinR)+cross(k_OinR,i_OinR_dot);
%     omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(v_MCR_t_dot,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
elseif strcmp(Flag_frame,'VVLH')
    j_OinR = -h_MCR_t/h_MCR_t_norm;
    k_OinR = -r_MCR_t/r_MCR_t_norm;
    i_OinR = cross(j_OinR,k_OinR);
    k_OinR_dot = - dot(v_MCR_t,i_OinR) / r_MCR_t_norm * i_OinR;
    j_OinR_dot = - dot(v_MCR_t_dot,j_OinR) * r_MCR_t_norm / h_MCR_t_norm * i_OinR;
    i_OinR_dot = r_MCR_t_norm/h_MCR_t_norm * dot(v_MCR_t_dot,j_OinR)*j_OinR + 1/r_MCR_t_norm * dot(v_MCR_t,i_OinR)*k_OinR;
%     omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(v_MCR_t_dot,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
else
    error('Wrong orbital frame flag')
end

% 坐标转换
TM_TCR2TCO = [i_OinR j_OinR k_OinR]';
r_TCO_c = TM_TCR2TCO*r_TCR_c;

if isempty(v_TCR_c)
    v_TCO_c = [];
else
    TM_TCR2TCO_dot = [i_OinR_dot,j_OinR_dot,k_OinR_dot]';
    v_TCO_c = TM_TCR2TCO*v_TCR_c + TM_TCR2TCO_dot*r_TCR_c;
%     v_TCO_c = TM_TCR2TCO*v_TCR_c + cross(omegaT,r_TCO_c);
end

end

