function F = maneuver_relmotion_eph(x0,dt,x0f_MCR_target,a0f_MCR_target,x0_TC_LVLH_chaser,rf_rel,aux)
x0_MCR_target = x0f_MCR_target(1,:);
xf_MCR_target = x0f_MCR_target(2,:);
a0_MCR_target = a0f_MCR_target(1,:);
af_MCR_target = a0f_MCR_target(2,:);
dv0 = x0(1:3);
x0_TC_LVLH_chaser(4:6) = x0_TC_LVLH_chaser(4:6)+dv0;

% chaser初值, LVLH 2 MCR
x0_MCR_chaser = T_TCO2TCR_eph(x0_TC_LVLH_chaser,x0_MCR_target,a0_MCR_target,'LVLH')+x0_MCR_target;

% 副星星历积分
xf_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,[0 dt],dt,aux);
% chaser, MCR 2 LVLH
xf_TC_rel = xf_MCR_chaser-xf_MCR_target;
xf_TC_LVLH = T_TCR2TCO_eph(xf_TC_rel,xf_MCR_target,af_MCR_target,'LVLH');

% 计算误差
F = xf_TC_LVLH(1:3)' - rf_rel(:);
