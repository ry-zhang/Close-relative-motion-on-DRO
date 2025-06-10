function [rvTC_VNC_c] = T_MCR2VNC(rvTC_MCR_c,rvMT_MCR_t)
%将跟踪航天器器相对位置rTC在月心（主星原点）旋转系MCR中的分量列阵转换到VNC系中
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH

%输入：
%    跟踪航天器相对位置在MCR中的分量列阵rTC_MCR
%    目标航天器位置速度rMT_MCR,vMT_MCR
%输出：
%    跟踪航天器相对位置在VNC中的分量列阵rTC_VNC

%轨道坐标系VNC
%    -原点位于目标航天器质心
%    -x轴沿着中心天体指向目标航天器
%    -z轴沿着目标航天器相对于中心天体的角速度方向
%月心旋转系MCR
%    -原点位于月心
%    -x轴沿月球指向地球方向
%    -z轴为瞬时月球相对于地球的角动量方向

row = size(rvTC_MCR_c,1);
rvTC_VNC_c = zeros(row,6);
for ii = 1:row
    rTC_MCR_c_ii = rvTC_MCR_c(ii,1:3); vTC_MCR_ii = rvTC_MCR_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_ii = rvMT_MCR_t(ii,4:6);
    [rTC_VNC_c_ii,vTC_VNC_c_ii] = MCR2VNC(rTC_MCR_c_ii',vTC_MCR_ii',rMT_MCR_t_ii',vMT_MCR_ii');
    rvTC_VNC_c(ii,:) = [rTC_VNC_c_ii',vTC_VNC_c_ii'];
end

end

function [rTC_VNC_c,vTC_VNC_c] = MCR2VNC(rTC_MCR_c,vTC_MCR_c,rMT_MCR_t,vMT_MCR_t)
hT = cross(rMT_MCR_t,vMT_MCR_t);
hT_norm = norm(hT);
rMT_norm = norm(rMT_MCR_t);
vMT_norm = norm(vMT_MCR_t);
omegaT = hT/rMT_norm^2;

%VNC坐标系定义
iVNC = vMT_MCR_t/vMT_norm;
jVNC = hT/hT_norm;
kVNC = cross(iVNC,jVNC);

%坐标转换
TM_VNC2MCR = [iVNC jVNC kVNC];
TM_MCR2VNC = TM_VNC2MCR';

rTC_VNC_c = TM_MCR2VNC*rTC_MCR_c;

vTC_VNC_c = TM_MCR2VNC*(vTC_MCR_c-cross(omegaT,rTC_MCR_c));
end

