function f = dist_DROformation(x0,xx_MCR_target,aa_MCR_target,tspan_sec,t_sample,r_MCRLVLH_rel_ref,scale_factor,aux)

% t_sample = linspace(tspan_sec(1),tspan_sec(2),500);
x0 = x0./scale_factor;
x0_TC_LVLH_chaser = [x0(1),x0(2),0,x0(3),x0(4),0];

x0_MCR_target = xx_MCR_target(1,:);
x0_MCR_chaser = T_TCO2TCR_eph(x0_TC_LVLH_chaser,x0_MCR_target,aa_MCR_target(1,:),'LVLH')+x0_MCR_target;

% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);

% 计算相对运动
rv_TC_MCR = xx_MCR_chaser-xx_MCR_target;
rv_TC_LVLH = T_TCR2TCO_eph(rv_TC_MCR,xx_MCR_target,aa_MCR_target,'LVLH');

% 距离参考轨道的最大距离
f = max(sqrt(sum((rv_TC_LVLH(:,[1,2,3])-r_MCRLVLH_rel_ref).^2,2))); % 三维最大距离,结果最好
% f = max(sqrt(sum((rv_TC_LVLH(:,[1,3])-r_MCRLVLH_rel_ref(:,[1,3])).^2,2))); % xy投影最大距离，结果次之
% f = mean(sqrt(sum((rv_TC_LVLH(:,[1,3])-r_MCRLVLH_rel_ref(:,[1,3])).^2,2)));% xy投影平均距离，结果最次
% disp([x0,f])