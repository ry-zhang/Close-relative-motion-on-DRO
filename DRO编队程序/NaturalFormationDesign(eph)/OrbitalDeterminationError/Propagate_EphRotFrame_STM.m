function [xx_MCR,a_MCR,STM_j2k] = Propagate_EphRotFrame_STM(x0_MCR,tspan_sec,t_sample,aux,UseParallel)
if nargin<5
    UseParallel = 0;
end

% 星历DRO初值 MCR → j2k
x0_j2k = T_Rot2ECJ2k(aux.jd0,x0_MCR,aux.C_Mat,'MCEMR',UseParallel);% km,km/s

% 星历积分设置
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
% 星历模型积分
sol = ode113(@(t, xx)eqom_geoMEMEJ2k_STM(t, xx, aux), ...
    tspan_sec, [x0_j2k,reshape(eye(6),1,36)] , options);

[xx_STM,a_STM] = deval(sol,t_sample);
xx_j2k = xx_STM(1:6,:)';
a_j2k = a_STM(1:6,:)';

STM_j2k = xx_STM(7:end,:)';

% DRO ECI → MCR
tt_jd = aux.jd0 + t_sample/86400;
[xx_MCR,a_MCR] = T_ECJ2k2Rot(tt_jd, xx_j2k, a_j2k(:,4:6),aux.C_Mat, 'MCEMR',UseParallel);


end