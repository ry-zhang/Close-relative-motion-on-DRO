function rvEme = BRC2Eme(jdate, rvRot, aux)
if size(rvRot,2)~=1
    rvRot = rvRot';
end
% mu_earth = aux.planet.mu(3);
% mu_moon = aux.planet.mu(10);
% mu_sun = aux.planet.mu(11);
% s -> day
% mu = mu_moon/(mu_earth+mu_moon);                                     %������
% EMRAT = 81.30056907419062; % DE430  % m_Earth/m_Moon
% mu = 1/(1+EMRAT);
mu = aux.mu;
rv_moon = ephEme_new(jdate, 10, 3 , aux.DE430);
r_moon = rv_moon(1:3);
v_moon = rv_moon(4:6);
a_moon = rv_moon(7:9);

 


drv =  dot(r_moon, v_moon);


crv = cross(r_moon, v_moon);
cra =  cross(r_moon, a_moon);

k = norm(r_moon);
k_dot = drv/k;

% k_2dot = (k*(dvv+dot(r_moon,a_moon)) - k_dot*(dot(r_moon,v_moon)))/k^2;
h = norm(crv);
h_dot = dot(crv,cra)/h;

% h_2dot =1/h^2* (h*((norm(cross(r_moon, a_moon)))^2+ dot(cross(r_moon, v_moon),  (cross(v_moon, a_moon)+cross(r_moon, j_moon)))) - h_dot*dot(cross(r_moon,v_moon),cross(r_moon,a_moon)));
%% ��ת����
e1 = r_moon/k;
e3 = crv/h;
e2 = cross(e3,e1);
C = [e1, e2, e3];
e1_dot = (k*v_moon - k_dot*r_moon)/k^2;
e3_dot = (h*cra-h_dot*crv)/h^2;
e2_dot = cross(e3_dot,e1) + cross(e3,e1_dot);

C_dot = [e1_dot, e2_dot, e3_dot];   
b = mu*r_moon;
b_dot = mu*v_moon;

% Omega_x = norm(r_moon)/(norm(h_moon))^2*(mu_sun/(norm(rMS))^3 -  mu_sun/(norm(r_sun))^3)*dot(r_sun, h_moon);
% Omega_y = 0;
% Omega_z = norm(h_moon)/(norm(r_moon))^2;
% Omega = [Omega_x; Omega_y ; Omega_z];                       %  ������תϵ��ʾ��Omega
% Omega = M_rot*Omega; 
% ������ת����(rot -> eme)
% b =zeros(3,1) ;
% b_dot =zeros(3,1);
omega = 1/aux.TU;
rEme =b+ k*C*rvRot(1:3);
vEme = b_dot+k_dot*C*rvRot(1:3) + k*C_dot*rvRot(1:3)+k*C*omega*rvRot(4:6) ;
 rvEme = [rEme; vEme];
end

