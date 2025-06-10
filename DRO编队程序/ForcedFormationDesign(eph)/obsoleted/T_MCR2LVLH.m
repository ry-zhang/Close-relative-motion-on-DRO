function [rvTC_LVLH_c] = T_MCR2LVLH(rvTC_MCR_c,rvMT_MCR_t)
%将跟踪航天器器相对位置rTC在月心（主星原点）旋转系MCR中的分量列阵转换到LVLH系中
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH
% LVLH坐标系即RIC

%输入：
%    跟踪航天器相对位置在MCR中的分量列阵rTC_MCR
%    目标航天器位置速度rMT_MCR,vMT_MCR
%输出：
%    跟踪航天器相对位置在LVLH中的分量列阵rTC_LVLH

%轨道坐标系LVLH
%    -原点位于目标航天器质心
%    -x轴沿着中心天体指向目标航天器
%    -z轴沿着目标航天器相对于中心天体的角速度方向
%月心旋转系MCR
%    -原点位于月心
%    -x轴沿月球指向地球方向
%    -z轴为瞬时月球相对于地球的角动量方向

row = size(rvTC_MCR_c,1);
rvTC_LVLH_c = zeros(row,6);
for ii = 1:row
    rTC_MCR_c_ii = rvTC_MCR_c(ii,1:3); vTC_MCR_ii = rvTC_MCR_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_ii = rvMT_MCR_t(ii,4:6);
    [rTC_LVLH_c_ii,vTC_LVLH_c_ii] = MCR2LVLH(rTC_MCR_c_ii',vTC_MCR_ii',rMT_MCR_t_ii',vMT_MCR_ii');
    rvTC_LVLH_c(ii,:) = [rTC_LVLH_c_ii',vTC_LVLH_c_ii'];
end

end

function [rTC_LVLH_c,vTC_LVLH_c] = MCR2LVLH(rTC_MCR_c,vTC_MCR_c,rMT_MCR_t,vMT_MCR_t)
hT = cross(rMT_MCR_t,vMT_MCR_t);
hT_norm = norm(hT);
rMT_norm = norm(rMT_MCR_t);
omegaT = hT/rMT_norm^2;

%LVLH坐标系定义
iLVLH = rMT_MCR_t/rMT_norm;
kLVLH = hT/hT_norm;
jLVLH = cross(kLVLH,iLVLH);

%坐标转换
TM_LVLH2MCR = [iLVLH jLVLH kLVLH];
TM_MCR2LVLH = TM_LVLH2MCR';

rTC_LVLH_c = TM_MCR2LVLH*rTC_MCR_c;

vTC_LVLH_c = TM_MCR2LVLH*(vTC_MCR_c-cross(omegaT,rTC_MCR_c));
end

