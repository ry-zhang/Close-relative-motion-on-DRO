function [xx_j2k,a_j2k] = Propagate_Ephj2k(x0_j2k,tspan_sec,t_sample,aux)

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

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
