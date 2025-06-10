function rvMC_MCR_c = T_VVLH2MCR(rvTC_VVLH_c,rvMT_MCR_t)
%轨道坐标系VVLH（非归一化）转化到月心（主星原点）旋转系MCR(归一化）
% v1 2020/8/20 fhl
% v2 2021/5/28 YCH
%输入：
%目标航天器位置速度rMT_MCR,vMT_MCR
%跟踪航天器位置速度rTC_VVLH,vTC_VVLH
%输出：
%跟踪航天器位置速度rTC_MCR,vTC_MCR
%轨道坐标系VVLH
%    -原点位于目标航天器质心
%    -y轴沿着目标航天器相对于中心天体的角速度的负方向
%    -z轴沿着目标航天器指向中心天体
%月心旋转系MCR
%    -原点位于月心
%    -x轴沿月球指向地球方向
%    -z轴为瞬时月球相对于地球的角动量方向
row = size(rvTC_VVLH_c,1);
rvMC_MCR_c = zeros(row,6);
for ii = 1:row
    rTC_VVLH_c_ii = rvTC_VVLH_c(ii,1:3); vTC_VVLH_c_ii = rvTC_VVLH_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_t_ii = rvMT_MCR_t(ii,4:6);
    [rTC_MCR_c_ii,vTC_MCR_c_ii] = VVLH2MCR(rTC_VVLH_c_ii',vTC_VVLH_c_ii',rMT_MCR_t_ii',vMT_MCR_t_ii');
    rvMC_MCR_c(ii,:) = [rTC_MCR_c_ii',vTC_MCR_c_ii'];
end

end

function [rMC_MCR_c,vMC_MCR_c] = VVLH2MCR(rTC_VVLH_c,vTC_VVLH_c,rMT_MCR_t,vMT_MCR_t)
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
rTC_MCR = TM_VVLH2MCR*rTC_VVLH_c;
rMC_MCR_c = rTC_MCR;
% rMC_MCR_c = rMT_MCR_t+rTC_MCR;

vTC_MCR = TM_VVLH2MCR*vTC_VVLH_c+cross(omegaT,rTC_MCR);
vMC_MCR_c = vTC_MCR;
% vMC_MCR_c = vMT_MCR_t+vTC_MCR;
end

