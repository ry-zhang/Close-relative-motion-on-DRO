function f = dist_DROformation_j2k(x0,xx_j2k_target,aa_j2k_target,tspan_sec,t_sample,r_ECILVLH_rel_ref,scale_factor,aux)

% t_sample = linspace(tspan_sec(1),tspan_sec(2),500);
x0 = x0./scale_factor;
x0_TC_j2kLVLH_chaser = [x0(1),x0(2),0,x0(3),x0(4),0];

x0_j2k_target = xx_j2k_target(1,:);
x0_j2k_chaser = T_TCO2TCR_eph(x0_TC_j2kLVLH_chaser,x0_j2k_target,aa_j2k_target(1,:),'LVLH')+x0_j2k_target;

% 副星星历积分
xx_j2k_chaser = Propagate_Ephj2k(x0_j2k_chaser,tspan_sec,t_sample,aux);

% 计算相对运动
rv_TC_j2k = xx_j2k_chaser-xx_j2k_target;
rv_TC_j2kLVLH = T_TCR2TCO_eph(rv_TC_j2k,xx_j2k_target,aa_j2k_target,'LVLH');

% 距离参考轨道的最大距离
f = max(sqrt(sum((rv_TC_j2kLVLH(:,[1,2,3])-r_ECILVLH_rel_ref).^2,2))); % 三维最大距离,结果最好
% f = max(sqrt(sum((rv_TC_LVLH(:,[1,3])-r_MCRLVLH_rel_ref(:,[1,3])).^2,2))); % xy投影最大距离，结果次之
% f = mean(sqrt(sum((rv_TC_LVLH(:,[1,3])-r_MCRLVLH_rel_ref(:,[1,3])).^2,2)));% xy投影平均距离，结果最次
% disp([x0,f])