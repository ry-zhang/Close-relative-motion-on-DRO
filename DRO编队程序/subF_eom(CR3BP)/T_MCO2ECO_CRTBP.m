function [xx_ECILVLH_rel,xx_ECI_target] = T_MCO2ECO_CRTBP(xx_MCRLVLH_rel,xx_MCR_target,t_sample,Flag_frame,mu)
% 月心旋转系(MCR)中的轨道坐标系转化到地心惯性系(ECI)中的轨道坐标系
% v1 2022/5/22 YCH

%输入：
%    跟踪航天器相对位置在MCR-轨道坐标系中的分量列阵xx_MCRLVLH_rel
%    目标航天器在MCR中的位置速度xx_MCR_target
%    坐标系标志，('VNC','LVLH','VVLH')
%输出：
%    跟踪航天器相对位置在TC中的分量列阵rv_TCI_c

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

xx_MCR_rel = T_TCO2TCR_CR3BP(xx_MCRLVLH_rel,xx_MCR_target,Flag_frame, mu);

% 将月球为中心天体的旋转坐标系转化到地球为中心天体的旋转坐标系 
xx_ECR_target = xx_MCR_target; xx_ECR_target(:,1:3) = xx_MCR_target(:,1:3)-[-1,0,0];
xx_ECR_rel = xx_MCR_rel;

% 地球为中心天体的惯性坐标系 下的轨道坐标系
xx_ECI_target = synodic2inertial(xx_ECR_target',t_sample)';
xx_ECI_rel = synodic2inertial(xx_ECR_rel',t_sample)';
xx_ECILVLH_rel = T_TCI2TCO_E_CR3BP(xx_ECI_rel,xx_ECI_target,t_sample,Flag_frame, mu);