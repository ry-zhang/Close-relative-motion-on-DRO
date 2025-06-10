function T_Mtx = get_Eme2RotJcoMtx(tao, aux)

jdate = aux.dep_t_jdtdb + tao*aux.TU/86400;
rv_12_eme = ephEme_new(jdate, 10 , 3 , aux.DE430); %% 月球相对于地球J2000惯性系的位置
r_12_eme = rv_12_eme(1 : 3);
v_12_eme = rv_12_eme(4 : 6);
a_12_eme = rv_12_eme(7:9);
omega_eme = cross(r_12_eme , v_12_eme);

h = norm(cross(r_12_eme , v_12_eme));
h_dot = dot(cross(r_12_eme,v_12_eme),cross(r_12_eme,a_12_eme))/h;

k = norm( r_12_eme);

k_dot = dot(r_12_eme, v_12_eme)/norm( r_12_eme);
% 计算旋转矩阵
e1 = r_12_eme / norm(r_12_eme);
e3 = omega_eme / norm(omega_eme);
e2 = cross(e3 , e1);
C = [e1 , e2 , e3];
e1_dot = (k*v_12_eme - k_dot* r_12_eme)/k^2;
e3_dot = (h*cross(r_12_eme , a_12_eme) - h_dot*cross(r_12_eme , v_12_eme))/h^2;
e2_dot = cross(e3_dot,e1) + cross(e3,e1_dot);
C_dot = [e1_dot,e2_dot,e3_dot];
mean_omega = 1/aux.TU;
T_Mtx = [k*C/aux.LU, zeros(3);
    (k_dot*C+k*C_dot)/aux.VU, k*C*mean_omega/aux.VU];
end

