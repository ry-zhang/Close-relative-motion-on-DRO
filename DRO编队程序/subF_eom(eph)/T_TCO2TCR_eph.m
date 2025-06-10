function rv_TCR_c = T_TCO2TCR_eph(rv_TCO_c,rv_MCR_t,a_MCR_t,Flag_frame)
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


% 判断输入的rv_TCO_c是否包含速度
[row,col] = size(rv_TCO_c);
if col == 3
    rv_TCR_c = zeros(row,3);
elseif col == 6
    rv_TCR_c = zeros(row,6);
else
    error('Wrong input size')
end

if isempty(a_MCR_t)
    a_MCR_t = zeros(row,3);
end

for ii = 1:row
    r_TCO_c_ii = rv_TCO_c(ii,1:3); 
    if col == 3
        v_TCO_c_ii = [];
    else
        v_TCO_c_ii = rv_TCO_c(ii,4:6);
    end
    r_MCR_t_ii = rv_MCR_t(ii,1:3); v_MCR_t_ii = rv_MCR_t(ii,4:6);a_MCR_t_ii = a_MCR_t(ii,:);
    [r_TCR_c_ii,v_TCR_c_ii] = TCO2TCR(r_TCO_c_ii',v_TCO_c_ii',r_MCR_t_ii',v_MCR_t_ii',a_MCR_t_ii',Flag_frame);
    rv_TCR_c(ii,:) = [r_TCR_c_ii',v_TCR_c_ii'];
end

end

function [r_TCR_c,v_TCR_c] = TCO2TCR(r_TCO_c,v_TCO_c,r_MCR_t,v_MCR_t,a_MCR_t,Flag_frame)

h_MCR_t = cross(r_MCR_t,v_MCR_t);
h_MCR_t_norm = norm(h_MCR_t);
v_MCR_t_norm = norm(v_MCR_t);
r_MCR_t_norm = norm(r_MCR_t);
% 坐标系定义
if strcmp(Flag_frame,'VNC')
    iVector = v_MCR_t/v_MCR_t_norm;
    jVector = h_MCR_t/h_MCR_t_norm;
    kVector = cross(iVector,jVector);
    omegaT = (dot(a_MCR_t,h_MCR_t)*cross(v_MCR_t,h_MCR_t)...
        + dot(a_MCR_t,h_MCR_t)*dot(r_MCR_t,v_MCR_t)*v_MCR_t...
        - dot(a_MCR_t,cross(v_MCR_t,h_MCR_t))*h_MCR_t)/h_MCR_t_norm^2/v_MCR_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    iVector = r_MCR_t/r_MCR_t_norm;
    kVector = h_MCR_t/h_MCR_t_norm;
    jVector = cross(kVector,iVector);
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
%     omegaT = h_MCR_t/r_MCR_t_norm^2;
elseif strcmp(Flag_frame,'VVLH')
    kVector = -r_MCR_t/r_MCR_t_norm;
    jVector = -h_MCR_t/h_MCR_t_norm;
    iVector = cross(jVector,kVector);
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
else
    error('Wrong orbital frame flag')
end

%坐标转换
TM_TCO2TCR = [iVector jVector kVector];
r_TCR_c = TM_TCO2TCR*r_TCO_c;
% r_TCR_c = r_TCR_c;
% rMC_TCR_c = r_TCR_t+rTCO_TCR;

if isempty(v_TCO_c)
    v_TCR_c = [];
else
    v_TCR_c = TM_TCO2TCR*v_TCO_c + cross(omegaT,r_TCR_c);
end
% v_TCR_c = v_TCR_c;
% vMC_TCR_c = v_TCR_t+vTCO_TCR;

end

