function delta_rf = maneuver_relmotion(x0,dt,x0_DRO,x0_REL,rf_rel,con)
dv0 = x0(1:3);
x0_REL(4:6) = x0_REL(4:6)+dv0;
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
[~,y] = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO, x0_REL], opts);
delta_rf = reshape(rf_rel(1:3),3,1) - y(end,7:9)';
% plot3(y(:,7)*para.r_norma,y(:,8)*para.r_norma,y(:,9)*para.r_norma);  axis equal
% xlabel('x [km]'); ylabel('y [km]')
delta_rf = delta_rf*10^2;



