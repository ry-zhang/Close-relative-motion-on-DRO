function [tf_sec,xf_MCEMR,tf_index] = find_thetatar(xx_MCEMR_AB,theta_tar)
% 依据一小时间隔的MCR下轨道 估算 积分到theta_tar附近所需的时间
% 如此，使每段优化时间段的起始端与末端位于theta_tar附近，以使轨道转移在观测范围内
theta_all = mod(atan2d(-xx_MCEMR_AB(:,3),xx_MCEMR_AB(:,2)),360); % 采用STK
theta_all_diff = theta_all-theta_tar;
nage_index = theta_all_diff(2:end).*theta_all_diff(1:end-1);
tf_index_all = find((nage_index<0)&(abs(theta_all_diff(2:end))<100)&abs(theta_all_diff(1:end-1))<100);
tf_index = tf_index_all(1);
tf_sec = xx_MCEMR_AB(tf_index,1);
xf_MCEMR = xx_MCEMR_AB(tf_index,2:end);
