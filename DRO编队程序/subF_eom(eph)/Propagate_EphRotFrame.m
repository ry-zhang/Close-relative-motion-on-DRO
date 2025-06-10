function [xx_MCR,a_MCR,aux_MCR_event] = Propagate_EphRotFrame(x0_MCR,tspan_sec,t_sample,aux,UseParallel,eventfunction,MaxStep)
if nargin<5
    UseParallel = 0;
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
elseif nargin == 5
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
elseif nargin == 6
    options = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',eventfunction);
elseif nargin == 7
    options = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',MaxStep,'Events',eventfunction);
end
% 星历DRO初值 MCR → j2k
x0_j2k = T_Rot2ECJ2k(aux.jd0+tspan_sec(1)/86400,x0_MCR,aux.C_Mat,'MCEMR',0);% km,km/s
tt_jd = aux.jd0 + t_sample/86400;

if tspan_sec(end) ~= tspan_sec(1)
    % 星历模型积分
    if isfield(aux,'jd0_2et')
        sol = ode113(@(t, xx)eom_SPICE_stm(t, xx, aux) , tspan_sec, x0_j2k , options);
    else
        sol = ode113(@(t, xx)eqom_geoMEMEJ2k(t, xx, aux), tspan_sec, x0_j2k , options);
    end
    [xx_j2k,a_j2k] = deval(sol,t_sample);
    xx_j2k = xx_j2k';
    a_j2k = a_j2k(4:6,:)';
else
    xx_j2k = x0_j2k;
    a_j2k = eomj2kMtx(xx_j2k,t_sample,aux);
end

% ECI → MCR
[xx_MCR,a_MCR] = T_ECJ2k2Rot(tt_jd, xx_j2k, a_j2k,aux.C_Mat, 'MCEMR',UseParallel);


if nargin >= 6
    [xe_j2k,ae_j2k] = deval(sol,sol.xe);
    xe_j2k = xe_j2k';
    ae_j2k = ae_j2k';
    % DRO ECI → MCR
    tt_jd = aux.jd0 + sol.xe/86400;
    aux_MCR_event.te = sol.xe;
    aux_MCR_event.ie = sol.ie;
    aux_MCR_event.xe_j2k = xe_j2k;
    aux_MCR_event.ae_j2k = ae_j2k;
    [aux_MCR_event.xe_MCR,aux_MCR_event.ae_MCR] = T_ECJ2k2Rot(tt_jd, xe_j2k, ae_j2k(:,4:6),aux.C_Mat, 'MCEMR',UseParallel);
else
    aux_MCR_event = [];
end

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