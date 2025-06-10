function [auxSFF,auxSTK] = sat_Coasting(auxSFF,auxSTK,alpha_tarAll)

theta_tarAll = alpha2theta(alpha_tarAll);
auxSFF.alpha_tarAll = alpha_tarAll;
auxSFF.theta_tarAll = theta_tarAll;

% 计算初始状态的theta与alpha
% 建立一个卫星auxSTK.sat_AB_coasting计算入轨之后coasting段的轨道
auxSTK.sat_AB_coasting = sat_initialize(auxSTK.scenario,auxSFF.x00_j2k_AB,...
    'CentralBody/Earth J2000',auxSTK.scenario.StartTime,20*86400,'sat_AB_coasting',auxSFF.PropagatorName);

t0_coast_sec = auxSTK.sat_AB_coasting.InitialState.OrbitEpoch;
tf = auxSTK.sat_AB_coasting.Propagate.StoppingConditions.Item('Duration').Properties.Trip;
xx_MCEMR_AB_coasting = PosVel_MCEMR(auxSTK.sat_AB_coasting.sat,t0_coast_sec,tf,3600);
theta_00 = mod(atan2d(-xx_MCEMR_AB_coasting(1,3),xx_MCEMR_AB_coasting(1,2)),360); % 采用STK

% 计算初始的alpha_0离哪个alpha_tar最近

auxSFF.theta_tar = theta_00+1e-7; % 强行设置theta与初始时刻相同，也即滑行段为0
auxSFF.alpha_tar = theta2alpha(auxSFF.theta_tar);

% theta_diff = mod(theta_tarAll-theta_00,360);
% [theta_sort,theta_sortindex] = sort(theta_diff);
% if theta_sort(1)>5 % 如果相位差大于5°
%     auxSFF.alpha_tar = alpha_tarAll(theta_sortindex(1));
%     auxSFF.theta_tar = theta_tarAll(1);
% else
%     auxSFF.alpha_tar = alpha_tarAll(theta_sortindex(2));
%     auxSFF.theta_tar = theta_tarAll(2);
% end

% 计算coasting段的递推时间并重置，计算coasting段在j2k下的轨道数据
[tf_coast_sec,xf_MCEMR_coast,tf_coast_index] = find_thetatar(xx_MCEMR_AB_coasting,auxSFF.theta_tar);
auxSTK.sat_AB_coasting.Propagate.StoppingConditions.Item('Duration').Properties.Trip = tf_coast_sec;
auxSTK.sat_AB_coasting.sat.Propagator.RunMCS;
xx_j2k_AB_coasting = PosVel_j2k(auxSTK.sat_AB_coasting.sat,t0_coast_sec,tf_coast_sec,3600);

% 保存coasting段的数据
auxSFF.coasting.t0 = 0;
auxSFF.coasting.tf = tf_coast_sec;
auxSFF.coasting.x0_MCEMR_AB = xf_MCEMR_coast;
auxSFF.coasting.x0_j2k_AB = xx_j2k_AB_coasting(end,2:end);
auxSFF.coasting.t_sample = xx_j2k_AB_coasting(:,1);
auxSFF.coasting.xx_MCEMR_AB = xx_MCEMR_AB_coasting(1:tf_coast_index,2:end);
auxSFF.coasting.xx_j2k_AB = xx_j2k_AB_coasting(:,2:end);
