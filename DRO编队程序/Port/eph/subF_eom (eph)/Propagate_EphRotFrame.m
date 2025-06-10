function [xx_MCR,a_MCR] = Propagate_EphRotFrame(x0_MCR,tspan_sec,t_sample,aux)
% 星历DRO初值 MCR → j2k
x0_j2k = T_Rot2ECJ2k(aux.jd0,x0_MCR,aux.C_Mat,'MCEMR');% km,km/s

% 星历积分设置
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% 星历模型积分
sol = ode113(@(t, xx)eqom_geoMEMEJ2k(t, xx, aux), tspan_sec, x0_j2k , options);

[xx_j2k,a_j2k] = deval(sol,t_sample);
% DRO ECI → MCR
tt_jd = aux.jd0 + t_sample/86400;
[xx_MCR,a_MCR] = T_ECJ2k2Rot(tt_jd, xx_j2k', a_j2k(4:6,:)',aux.C_Mat, 'MCEMR');

% var_tested = xx_MCR;
% % var_tested = [xx_MCR(:,4:6),a_MCR];
% R = diff(var_tested(2:end,1:3))./var_tested(2:end-1,4:6)/diff(t_sample(1:2));
% plot(t_sample(3:end)'./86400,R,'.','LineWidth',1.5); 
% xlabel('t [day]')
% xlim([0,max(t_sample)./86400]); 
% ylim([0.95,1.05])
% legend('R_x','R_y','R_z')
% set(gca,'FontSize',15,'fontname','times new roman');
end