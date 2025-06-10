function Phi_Rot2ECJ2k = T_TCO2TCR_eph_phi(rv_MCR_t,a_MCR_t,Flag_frame)
% �������ϵ(target-centered orbital,TCO)ת��������ԭ����תϵ(target-centered rotational,TCR)
% v1 2020/8/20 fhl
% v2 2021/6/6 YCH ������ٶȣ����������ϵ��ת��������һ��

%���룺
%    ���ٺ��������λ����TCR�еķ�������rv_TCO_c
%    Ŀ�꺽����λ���ٶ�rvMT_MCR_t
%    ����ϵ��־��('VNC','LVLH','VVLH')
%�����
%    ���ٺ��������λ����VNC�еķ�������rTCO_c

%�������ϵVNC
%    -ԭ��λ��Ŀ�꺽��������
%    -x������Ŀ�꺽�������ٶȷ���
%    -y������Ŀ�꺽�����������������Ľ��ٶȷ��� 
%�������ϵLVLH
%    -ԭ��λ��Ŀ�꺽��������
%    -x��������������ָ��Ŀ�꺽����
%    -z������Ŀ�꺽�����������������Ľ��ٶȷ���
%�������ϵVVLH
%    -ԭ��λ��Ŀ�꺽��������
%    -y������Ŀ�꺽�����������������Ľ��ٶȵĸ�����
%    -z������Ŀ�꺽����ָ����������

% �ؼ�����rv_MCR_t,a_MCR_t������λ��ͬһ�������ᶨ�������ϵ��
% rv_TCR_c���ڵ�Ŀ������ϵ���Ǹ��������r_MCR_t,v_MCR_tʸ�������;
% ��˿�����������Ŀ������ϵ��������תϵ��j2k�ȵȣ�

rv_MCR_t = rv_MCR_t(:);
a_MCR_t = a_MCR_t(:);

r_MCR_t = rv_MCR_t(1:3); v_MCR_t = rv_MCR_t(4:6);

h_MCR_t = cross(r_MCR_t,v_MCR_t);
h_MCR_t_norm = norm(h_MCR_t);
v_MCR_t_norm = norm(v_MCR_t);
r_MCR_t_norm = norm(r_MCR_t);

% ����ϵ����
if strcmp(Flag_frame,'VNC')
    iVector = v_MCR_t/v_MCR_t_norm;
    jVector = h_MCR_t/h_MCR_t_norm;
    kVector = cross(iVector,jVector);
    iVector_dot = - dot(v_MCR_t,a_MCR_t) / v_MCR_t_norm^3 * v_MCR_t + a_MCR_t/v_MCR_t_norm;
    jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * cross(jVector,r_MCR_t/r_MCR_t_norm);
    kVector_dot = cross(iVector_dot,jVector)+cross(iVector,jVector_dot);
%     omegaT = (dot(a_MCR_t,h_MCR_t)*cross(v_MCR_t,h_MCR_t)...
%         + dot(a_MCR_t,h_MCR_t)*dot(r_MCR_t,v_MCR_t)*v_MCR_t...
%         - dot(a_MCR_t,cross(v_MCR_t,h_MCR_t))*h_MCR_t)/h_MCR_t_norm^2/v_MCR_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    iVector = r_MCR_t/r_MCR_t_norm;
    kVector = h_MCR_t/h_MCR_t_norm;
    jVector = cross(kVector,iVector);
    iVector_dot = dot(v_MCR_t,jVector) / r_MCR_t_norm * jVector;
    kVector_dot = - dot(a_MCR_t,kVector) * r_MCR_t_norm / h_MCR_t_norm * jVector;
    jVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,kVector)*kVector - 1/r_MCR_t_norm * dot(v_MCR_t,jVector)*iVector;
%     omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
elseif strcmp(Flag_frame,'VVLH')
    kVector = -r_MCR_t/r_MCR_t_norm;
    jVector = -h_MCR_t/h_MCR_t_norm;
    iVector = cross(jVector,kVector);
    kVector_dot = - dot(v_MCR_t,iVector) / r_MCR_t_norm * iVector;
    jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * iVector;
    iVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,jVector)*jVector + 1/r_MCR_t_norm * dot(v_MCR_t,iVector)*kVector;
%     omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
else
    error('Wrong orbital frame flag')
end

%����ת��
TM_TCO2TCR = [iVector jVector kVector];
TM_TCO2TCR_dot = [iVector_dot jVector_dot kVector_dot];
Phi_Rot2ECJ2k = [TM_TCO2TCR, zeros(3); TM_TCO2TCR_dot, TM_TCO2TCR];
end

