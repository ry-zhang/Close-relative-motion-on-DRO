function Phiv_det = xysingularTM(dt,t0,x0_DRO,con)
% 状态转移矩阵求xy方向转移奇点

opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

tf = t0+dt;
[~,sol_DRO_rel] = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[t0 tf], [x0_DRO', zeros(1,6),reshape(eye(6,6),1,36)], opts);
Phi = reshape(sol_DRO_rel(end,13:end),6,6);
Phiv_det = det(Phi(1:2,4:5));


