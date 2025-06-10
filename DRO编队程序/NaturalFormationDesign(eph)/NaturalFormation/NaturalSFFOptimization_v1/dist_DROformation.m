function f = dist_DROformation(x0,xx_MCR_target,a_MCR_target,tspan_sec,t_sample,rRel_ref,scale_factor,aux)

% t_sample = linspace(tspan_sec(1),tspan_sec(2),500);
x0 = x0./scale_factor;
x0_TC_VVLH_chaser = [x0(1),0,x0(2),x0(3),0,x0(4)];

x0_MCR_target = xx_MCR_target(1,:);
x0_MCR_chaser = T_TCO2TCR_eph(x0_TC_VVLH_chaser,x0_MCR_target,a_MCR_target(1,:),'VVLH')+x0_MCR_target;

% 副星星历积分
xx_MCR_chaser = Propagate_EphRotFrame(x0_MCR_chaser,tspan_sec,t_sample,aux);

% 计算相对运动
rvTC_MCR = xx_MCR_chaser-xx_MCR_target;
rvTC_VVLH = T_TCR2TCO_eph(rvTC_MCR,xx_MCR_target,a_MCR_target,'VVLH');

% 距离参考轨道的最大距离
rRel_ref = [rRel_ref(:,1),zeros(size(rRel_ref(:,1))),rRel_ref(:,2)];
f = max(sqrt(sum((rvTC_VVLH(:,[1,2,3])-rRel_ref).^2,2))); % 三维
% f = max(sqrt(sum((rvTC_VVLH(:,[1,3])-rRel_ref).^2,2))); %只看xy投影
% f = mean(sqrt(sum((rvTC_VVLH(:,[1,3])-rRel_ref).^2,2)));
% disp([x0,f])