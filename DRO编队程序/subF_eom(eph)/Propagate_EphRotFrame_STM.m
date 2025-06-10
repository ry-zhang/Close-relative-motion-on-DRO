function [xx_MCR,a_MCR,phi_MCR,DROtheta] = Propagate_EphRotFrame_STM(x0_MCR,tspan_sec,t_sample,aux,UseParallel,eventfunction)
if nargin<5
    UseParallel = 0;
    % 星历积分设置
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
elseif nargin == 6
    options = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',eventfunction);
end
% 星历DRO初值 MCR → j2k
x0_j2k = T_Rot2ECJ2k(aux.jd0,x0_MCR,aux.C_Mat,'MCEMR',UseParallel);% km,km/s


% 星历模型积分
phi0 = eye(6);
if isfield(aux,'jd0_2et')
    sol = ode113(@(t, xx)eom_SPICE_stm(t, xx, aux) , tspan_sec, [x0_j2k, phi0(:)'] , options);
else
    sol = ode113(@(t, xx)eqom_geoMEMEJ2k_STM(t, xx, aux), tspan_sec, [x0_j2k, phi0(:)'] , options);
end

[xx_j2k,a_j2k] = deval(sol,t_sample);
phi_j2k = reshape(xx_j2k(7:end,end),6,6);
xx_j2k = xx_j2k(1:6,:)';
a_j2k = a_j2k(1:6,:)';


% DRO ECI → MCR
tt_jd = aux.jd0 + t_sample/86400;
[xx_MCR,a_MCR] = T_ECJ2k2Rot(tt_jd, xx_j2k, a_j2k(:,4:6),aux.C_Mat, 'MCEMR',UseParallel);

% STM j2k -> MCR
Phi_Rot2ECJ2k_0 = T_Rot2ECJ2k_phi(aux.jd0, aux.C_Mat, 'MCEMR');
Phi_Rot2ECJ2k_f = T_Rot2ECJ2k_phi(tt_jd(end), aux.C_Mat, 'MCEMR');
phi_MCR = Phi_Rot2ECJ2k_f^(-1)*phi_j2k*Phi_Rot2ECJ2k_0;

if nargin == 6
    [xe_j2k,ae_j2k] = deval(sol,sol.xe);
    xe_j2k = xe_j2k';
    ae_j2k = ae_j2k';
    % DRO ECI → MCR
    tt_jd = aux.jd0 + sol.xe/86400;
    DROtheta.te = sol.xe;
    DROtheta.ie = sol.ie;
    [DROtheta.xe_MCR,DROtheta.ae_MCR] = T_ECJ2k2Rot(tt_jd, xe_j2k, ae_j2k(:,4:6),aux.C_Mat, 'MCEMR',UseParallel);
else
    DROtheta = [];
end

end