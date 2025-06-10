
%% j2k STM 数值验证
phi_num = zeros(6);
dx = 1e-4;
for ii = 1:6
    x0_j2k_d = x0_j2k; x0_j2k_d(ii) = x0_j2k(ii) + dx;
    sol_d = ode113(@(t, xx)eqom_geoMEMEJ2k_STM(t, xx, aux), tspan_sec, [x0_j2k_d, phi0(:)'], options);
    phi_num(:,ii) = (sol_d.y(1:6,end)-sol.y(1:6,end))/dx;
end

%% MCR STM 数值验证
phi_num = zeros(6);
dx = 1e-4;
for ii = 1:6
    x0_MCR_chaser_d = x0_MCR_chaser; x0_MCR_chaser_d(ii) = x0_MCR_chaser(ii) + dx;
    xf_MCR_chaser_d = Propagate_EphRotFrame_STM(x0_MCR_chaser_d,[0 dt],dt,aux);
    phi_num(:,ii) = (xf_MCR_chaser_d-xf_MCR_chaser_d)/dx;
end

%% TCO STM 数值验证
phi_num = zeros(6);
dx = 1e-5;
for ii = 1:6
    x0_TC_VVLH_chaser_d = x0_TC_VVLH_chaser; 
    x0_TC_VVLH_chaser_d(ii) = x0_TC_VVLH_chaser(ii) + dx;
    % chaser初值, VVLH 2 MCR
    x0_MCR_chaser_d = T_TCO2TCR_eph(x0_TC_VVLH_chaser_d,x0_MCR_target,a0_MCR_target,'VVLH')+x0_MCR_target;

    % 副星星历积分
    xf_MCR_chaser_d = Propagate_EphRotFrame_STM(x0_MCR_chaser_d,[0 dt],dt,aux);


    % chaser, MCR 2 VVLH
    xf_TC_rel_d = xf_MCR_chaser_d-xf_MCR_target;
    xf_TC_VVLH_d = T_TCR2TCO_eph(xf_TC_rel_d,xf_MCR_target,af_MCR_target,'VVLH');
    
    phi_num(:,ii) = (xf_TC_VVLH_d-xf_TC_VVLH)/dx;
end

%% maneuver_relmotion_eph_STM 数值验证
phi_num = zeros(3);
[F,J] = maneuver_relmotion_eph_STM(x0,dt,x0f_MCR_target,a0f_MCR_target,x0_REL_loop,rf_REL,aux);
dx = 1e-6;
for ii = 1:3
    x0_d = x0; x0_d(ii) = x0(ii) + dx;
    Fd = maneuver_relmotion_eph_STM(x0_d,dt,x0f_MCR_target,a0f_MCR_target,x0_REL_loop,rf_REL,aux);
    
    phi_num(:,ii) = (Fd-F)/dx;
end