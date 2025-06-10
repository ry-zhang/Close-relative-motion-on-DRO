function Phi_Rot2ECJ2k = T_TCO2TCR_eph_phi(rv_MCR_t,a_MCR_t,Flag_frame)
% 轨道坐标系(target-centered orbital,TCO)转化到主星原点旋转系(target-centered rotational,TCR)
% v1 2020/8/20 fhl
% v2 2021/6/6 YCH 加入角速度，将多个坐标系的转换集成在一起

%输入：
%    跟踪航天器相对位置在TCR中的分量列阵rv_TCO_c
%    目标航天器位置速度rvMT_MCR_t
%    坐标系标志，('VNC','LVLH','VVLH')
%输出：
%    跟踪航天器相对位置在VNC中的分量列阵rTCO_c

%轨道坐标系VNC
%    -原点位于目标航天器质心
%    -x轴沿着目标航天器的速度方向
%    -y轴沿着目标航天器相对于中心天体的角速度方向 
%轨道坐标系LVLH
%    -原点位于目标航天器质心
%    -x轴沿着中心天体指向目标航天器
%    -z轴沿着目标航天器相对于中心天体的角速度方向
%轨道坐标系VVLH
%    -原点位于目标航天器质心
%    -y轴沿着目标航天器相对于中心天体的角速度的负方向
%    -z轴沿着目标航天器指向中心天体

% 关键在于rv_MCR_t,a_MCR_t两者需位于同一个坐标轴定义的坐标系下
% rv_TCR_c所在的目标坐标系即是根据输入的r_MCR_t,v_MCR_t矢量定义的;
% 因此可适用于任意目标坐标系（月心旋转系，j2k等等）

rv_MCR_t = rv_MCR_t(:);
a_MCR_t = a_MCR_t(:);

r_MCR_t = rv_MCR_t(1:3); v_MCR_t = rv_MCR_t(4:6);

h_MCR_t = cross(r_MCR_t,v_MCR_t);
h_MCR_t_norm = norm(h_MCR_t);
v_MCR_t_norm = norm(v_MCR_t);
r_MCR_t_norm = norm(r_MCR_t);

% 坐标系定义
if strcmp(Flag_frame,'VNC')
    iVector = v_MCR_t/v_MCR_t_norm;
    jVector = h_MCR_t/h_MCR_t_norm;
    kVector = cross(iVector,jVector);
    iVector_dot = - dot(v_MCR_t,a_MCR_t) / v_MCR_t_norm^3 * v_MCR_t + a_MCR_t/v_MCR_t_norm;
    jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * cross(jVector,r_MCR_t/r_MCR_t_norm);
    kVector_dot = cross(iVector_dot,jVector)+cross(iVector,jVector_dot);
%     omegaT = (dot(a_MCR_t,h_MCR_t)*cross(v_MCR_t,h_MCR_t)...
%         + dot(a_MCR_t,h_MCR_t)*dot(r_MCR_t,v_MCR_t)*v_MCR_t...
%         - dot(a_MCR_t,cross(v_MCR_t,h_MCR_t))*h_MCR_t)/h_MCR_t_norm^2/v_MCR_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    iVector = r_MCR_t/r_MCR_t_norm;
    kVector = h_MCR_t/h_MCR_t_norm;
    jVector = cross(kVector,iVector);
    iVector_dot = dot(v_MCR_t,jVector) / r_MCR_t_norm * jVector;
    kVector_dot = - dot(a_MCR_t,kVector) * r_MCR_t_norm / h_MCR_t_norm * jVector;
    jVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,kVector)*kVector - 1/r_MCR_t_norm * dot(v_MCR_t,jVector)*iVector;
%     omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
elseif strcmp(Flag_frame,'VVLH')
    kVector = -r_MCR_t/r_MCR_t_norm;
    jVector = -h_MCR_t/h_MCR_t_norm;
    iVector = cross(jVector,kVector);
    kVector_dot = - dot(v_MCR_t,iVector) / r_MCR_t_norm * iVector;
    jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * iVector;
    iVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,jVector)*jVector + 1/r_MCR_t_norm * dot(v_MCR_t,iVector)*kVector;
%     omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
else
    error('Wrong orbital frame flag')
end

%坐标转换
TM_TCO2TCR = [iVector jVector kVector];
TM_TCO2TCR_dot = [iVector_dot jVector_dot kVector_dot];
Phi_Rot2ECJ2k = [TM_TCO2TCR, zeros(3); TM_TCO2TCR_dot, TM_TCO2TCR];
end

