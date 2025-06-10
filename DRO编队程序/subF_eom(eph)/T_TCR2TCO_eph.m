function [rv_TCO_c] = T_TCR2TCO_eph(rv_TCR_c,rv_MCR_t,a_MCR_t,Flag_frame)
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

% 关键在于rv_TCR_c,rv_MCR_t,a_MCR_t三者需位于同一个坐标轴定义的坐标系下
% 原坐标系即是根据输入的r_MCR_t,v_MCR_t矢量定义的;
% 因此可适用于任意原坐标系（月心旋转系，j2k等等）
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

if isempty(a_MCR_t)
    a_MCR_t = zeros(row,3);
end

for ii = 1:row
    r_TCR_c_ii = rv_TCR_c(ii,1:3); 
    if col == 3
        v_TCR_c_ii = [];
    else
        v_TCR_c_ii = rv_TCR_c(ii,4:6);
    end
    r_MCR_t_ii = rv_MCR_t(ii,1:3); v_MCR_ii = rv_MCR_t(ii,4:6);a_MCR_t_ii = a_MCR_t(ii,:);
    [r_TC_c_ii,v_TC_c_ii] = TCR2TCO(r_TCR_c_ii',v_TCR_c_ii',r_MCR_t_ii',v_MCR_ii',a_MCR_t_ii',Flag_frame);
    rv_TCO_c(ii,:) = [r_TC_c_ii',v_TC_c_ii'];
end

end

function [r_TCO_c,v_TCO_c] = TCR2TCO(r_TCR_c,v_TCR_c,r_MCR_t,v_MCR_t,a_MCR_t,Flag_frame)

r_MCR_t_norm = norm(r_MCR_t);
h_MCR_t = cross(r_MCR_t,v_MCR_t);
h_MCR_t_norm = norm(h_MCR_t);
v_MCR_t_norm = norm(v_MCR_t);

% 坐标系定义
if strcmp(Flag_frame,'VNC')
    iVector = v_MCR_t/v_MCR_t_norm;
    jVector = h_MCR_t/h_MCR_t_norm;
    kVector = cross(iVector,jVector);
%     iVector_dot = - dot(v_MCR_t,a_MCR_t) / v_MCR_t_norm^3 * v_MCR_t + a_MCR_t/v_MCR_t_norm;
%     jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * cross(jVector,r_MCR_t/r_MCR_t_norm);
%     kVector_dot = cross(iVector_dot,jVector)+cross(iVector,jVector_dot);
    omegaT = (dot(a_MCR_t,h_MCR_t)*cross(v_MCR_t,h_MCR_t)...
        + dot(a_MCR_t,h_MCR_t)*dot(r_MCR_t,v_MCR_t)*v_MCR_t...
        - dot(a_MCR_t,cross(v_MCR_t,h_MCR_t))*h_MCR_t)/h_MCR_t_norm^2/v_MCR_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    iVector = r_MCR_t/r_MCR_t_norm;
    kVector = h_MCR_t/h_MCR_t_norm;
    jVector = cross(kVector,iVector);
%     iVector_dot = dot(v_MCR_t,jVector) / r_MCR_t_norm * jVector;
%     kVector_dot = - dot(a_MCR_t,kVector) * r_MCR_t_norm / h_MCR_t_norm * jVector;
%     jVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,kVector)*kVector - 1/r_MCR_t_norm * dot(v_MCR_t,jVector)*iVector;
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
elseif strcmp(Flag_frame,'VVLH')
    kVector = -r_MCR_t/r_MCR_t_norm;
    jVector = -h_MCR_t/h_MCR_t_norm;
    iVector = cross(jVector,kVector);
%     kVector_dot = - dot(v_MCR_t,iVector) / r_MCR_t_norm * iVector;
%     jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * iVector;
%     iVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,jVector)*jVector + 1/r_MCR_t_norm * dot(v_MCR_t,iVector)*kVector;
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
else
    error('Wrong orbital frame flag')
end

%坐标转换
TM_TCR2TCO = [iVector jVector kVector]';
r_TCO_c = TM_TCR2TCO*r_TCR_c;
% r_TCR_c = r_TCR_c;
% rMC_TCR_c = r_TCR_t+rTCO_TCR;

if isempty(v_TCR_c)
    v_TCO_c = [];
else
    v_TCO_c = TM_TCR2TCO*v_TCR_c - cross(TM_TCR2TCO*omegaT,r_TCO_c);
%     TM_TCR2TCO_dot = [iVector_dot,jVector_dot,kVector_dot]';
%     v_TCO_c = TM_TCR2TCO*v_TCR_c + TM_TCR2TCO_dot*r_TCR_c;
end
% v_TCR_c = v_TCR_c;
% vMC_TCR_c = v_TCR_t+vTCO_TCR;

end

