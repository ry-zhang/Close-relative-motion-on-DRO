function [rvTC_VVLH_c] = T_MCR2VVLH(rvTC_MCR_c,rvMT_MCR_t)
%将跟踪航天器器相对位置rTC在月心（主星原点）旋转系MCR中的分量列阵转换到VVLH系中
% v1 2020/8/20 fhl
% v2 2021/5/28 YCH
%输入：
%    跟踪航天器相对位置在MCR中的分量列阵rTC_MCR
%    目标航天器位置速度rMT_MCR,vMT_MCR
%输出：
%    跟踪航天器相对位置在VVLH中的分量列阵rTC_VVLH
%轨道坐标系VVLH
%    -原点位于目标航天器质心
%    -y轴沿着目标航天器相对于中心天体的角速度的负方向
%    -z轴沿着目标航天器指向中心天体
%月心旋转系MCR
%    -原点位于月心
%    -x轴沿月球指向地球方向
%    -z轴为瞬时月球相对于地球的角动量方向

row = size(rvTC_MCR_c,1);
rvTC_VVLH_c = zeros(row,6);
for ii = 1:row
    rTC_MCR_c_ii = rvTC_MCR_c(ii,1:3); vTC_MCR_ii = rvTC_MCR_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_ii = rvMT_MCR_t(ii,4:6);
    [rTC_VVLH_c_ii,vTC_VVLH_c_ii] = MCR2VVLH(rTC_MCR_c_ii',vTC_MCR_ii',rMT_MCR_t_ii',vMT_MCR_ii');
    rvTC_VVLH_c(ii,:) = [rTC_VVLH_c_ii',vTC_VVLH_c_ii'];
end

end

function [rTC_VVLH_c,vTC_VVLH_c] = MCR2VVLH(rTC_MCR_c,vTC_MCR_c,rMT_MCR_t,vMT_MCR_t)
hT = cross(rMT_MCR_t,vMT_MCR_t);
hT_norm = norm(hT);
rMT_norm = norm(rMT_MCR_t);
omegaT = hT/rMT_norm^2;

%VVLH坐标系定义
kVVLH = -rMT_MCR_t/rMT_norm;
jVVLH = -hT/hT_norm;
iVVLH = cross(jVVLH,kVVLH);

%坐标转换
TM_VVLH2MCR = [iVVLH jVVLH kVVLH];
TM_MCR2VVLH = TM_VVLH2MCR';

rTC_VVLH_c = TM_MCR2VVLH*rTC_MCR_c;

vTC_VVLH_c = TM_MCR2VVLH*(vTC_MCR_c-cross(omegaT,rTC_MCR_c));
end

