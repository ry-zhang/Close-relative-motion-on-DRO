function [rv_TCO_c] = T_TCR2TCO_eph(rv_TCR_c,rv_MCR_t,a_MCR_t,Flag_frame)
% ����ԭ����תϵ(target-centered rotational,TCR)ת�����������ϵ(target-centered orbital,TCO)
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH

%���룺
%    ���ٺ��������λ����MCR�еķ�������rv_MCR_t
%    Ŀ�꺽����λ���ٶ�rMT_MCR,vMT_MCR
%    ����ϵ��־��('VNC','LVLH','VVLH')
%�����
%    ���ٺ��������λ����TC�еķ�������rv_TCR_c

%�������ϵVNC
%    -ԭ��λ��Ŀ�꺽��������
%    -x������Ŀ�꺽�������ٶȷ���
%    -y������Ŀ�꺽�����������������Ľ��ٶȷ��� 
%�������ϵLVLH(��RIC)
%    -ԭ��λ��Ŀ�꺽��������
%    -x��������������ָ��Ŀ�꺽����
%    -z������Ŀ�꺽�����������������Ľ��ٶȷ���
%�������ϵVVLH
%    -ԭ��λ��Ŀ�꺽��������
%    -y������Ŀ�꺽�����������������Ľ��ٶȵĸ�����
%    -z������Ŀ�꺽����ָ����������

% �ؼ�����rv_TCR_c,rv_MCR_t,a_MCR_t������λ��ͬһ�������ᶨ�������ϵ��
% ԭ����ϵ���Ǹ��������r_MCR_t,v_MCR_tʸ�������;
% ��˿�����������ԭ����ϵ��������תϵ��j2k�ȵȣ�
%    -ԭ��λ������
%    -x��������ָ�������
%    -z��Ϊ˲ʱ��������ڵ���ĽǶ�������

% �ж������rv_TCR_c�Ƿ�����ٶ�
[row,col] = size(rv_TCR_c);
if col == 3
    rv_TCO_c = zeros(row,3);
elseif col == 6
    rv_TCO_c = zeros(row,6);
else
    error('Wrong input size')
end

if isempty(a_MCR_t)
    a_MCR_t = zeros(row,3);
end

for ii = 1:row
    r_TCR_c_ii = rv_TCR_c(ii,1:3); 
    if col == 3
        v_TCR_c_ii = [];
    else
        v_TCR_c_ii = rv_TCR_c(ii,4:6);
    end
    r_MCR_t_ii = rv_MCR_t(ii,1:3); v_MCR_ii = rv_MCR_t(ii,4:6);a_MCR_t_ii = a_MCR_t(ii,:);
    [r_TC_c_ii,v_TC_c_ii] = TCR2TCO(r_TCR_c_ii',v_TCR_c_ii',r_MCR_t_ii',v_MCR_ii',a_MCR_t_ii',Flag_frame);
    rv_TCO_c(ii,:) = [r_TC_c_ii',v_TC_c_ii'];
end

end

function [r_TCO_c,v_TCO_c] = TCR2TCO(r_TCR_c,v_TCR_c,r_MCR_t,v_MCR_t,a_MCR_t,Flag_frame)

r_MCR_t_norm = norm(r_MCR_t);
h_MCR_t = cross(r_MCR_t,v_MCR_t);
h_MCR_t_norm = norm(h_MCR_t);
v_MCR_t_norm = norm(v_MCR_t);

% ����ϵ����
if strcmp(Flag_frame,'VNC')
    iVector = v_MCR_t/v_MCR_t_norm;
    jVector = h_MCR_t/h_MCR_t_norm;
    kVector = cross(iVector,jVector);
%     iVector_dot = - dot(v_MCR_t,a_MCR_t) / v_MCR_t_norm^3 * v_MCR_t + a_MCR_t/v_MCR_t_norm;
%     jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * cross(jVector,r_MCR_t/r_MCR_t_norm);
%     kVector_dot = cross(iVector_dot,jVector)+cross(iVector,jVector_dot);
    omegaT = (dot(a_MCR_t,h_MCR_t)*cross(v_MCR_t,h_MCR_t)...
        + dot(a_MCR_t,h_MCR_t)*dot(r_MCR_t,v_MCR_t)*v_MCR_t...
        - dot(a_MCR_t,cross(v_MCR_t,h_MCR_t))*h_MCR_t)/h_MCR_t_norm^2/v_MCR_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    iVector = r_MCR_t/r_MCR_t_norm;
    kVector = h_MCR_t/h_MCR_t_norm;
    jVector = cross(kVector,iVector);
%     iVector_dot = dot(v_MCR_t,jVector) / r_MCR_t_norm * jVector;
%     kVector_dot = - dot(a_MCR_t,kVector) * r_MCR_t_norm / h_MCR_t_norm * jVector;
%     jVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,kVector)*kVector - 1/r_MCR_t_norm * dot(v_MCR_t,jVector)*iVector;
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
elseif strcmp(Flag_frame,'VVLH')
    kVector = -r_MCR_t/r_MCR_t_norm;
    jVector = -h_MCR_t/h_MCR_t_norm;
    iVector = cross(jVector,kVector);
%     kVector_dot = - dot(v_MCR_t,iVector) / r_MCR_t_norm * iVector;
%     jVector_dot = - dot(a_MCR_t,jVector) * r_MCR_t_norm / h_MCR_t_norm * iVector;
%     iVector_dot = r_MCR_t_norm/h_MCR_t_norm * dot(a_MCR_t,jVector)*jVector + 1/r_MCR_t_norm * dot(v_MCR_t,iVector)*kVector;
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
else
    error('Wrong orbital frame flag')
end

%����ת��
TM_TCR2TCO = [iVector jVector kVector]';
r_TCO_c = TM_TCR2TCO*r_TCR_c;
% r_TCR_c = r_TCR_c;
% rMC_TCR_c = r_TCR_t+rTCO_TCR;

if isempty(v_TCR_c)
    v_TCO_c = [];
else
    v_TCO_c = TM_TCR2TCO*v_TCR_c - cross(TM_TCR2TCO*omegaT,r_TCO_c);
%     TM_TCR2TCO_dot = [iVector_dot,jVector_dot,kVector_dot]';
%     v_TCO_c = TM_TCR2TCO*v_TCR_c + TM_TCR2TCO_dot*r_TCR_c;
end
% v_TCR_c = v_TCR_c;
% vMC_TCR_c = v_TCR_t+vTCO_TCR;

end

