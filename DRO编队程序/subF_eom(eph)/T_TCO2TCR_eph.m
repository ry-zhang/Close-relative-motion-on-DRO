function rv_TCR_c = T_TCO2TCR_eph(rv_TCO_c,rv_MCR_t,a_MCR_t,Flag_frame)
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


% �ж������rv_TCO_c�Ƿ�����ٶ�
[row,col] = size(rv_TCO_c);
if col == 3
    rv_TCR_c = zeros(row,3);
elseif col == 6
    rv_TCR_c = zeros(row,6);
else
    error('Wrong input size')
end

if isempty(a_MCR_t)
    a_MCR_t = zeros(row,3);
end

for ii = 1:row
    r_TCO_c_ii = rv_TCO_c(ii,1:3); 
    if col == 3
        v_TCO_c_ii = [];
    else
        v_TCO_c_ii = rv_TCO_c(ii,4:6);
    end
    r_MCR_t_ii = rv_MCR_t(ii,1:3); v_MCR_t_ii = rv_MCR_t(ii,4:6);a_MCR_t_ii = a_MCR_t(ii,:);
    [r_TCR_c_ii,v_TCR_c_ii] = TCO2TCR(r_TCO_c_ii',v_TCO_c_ii',r_MCR_t_ii',v_MCR_t_ii',a_MCR_t_ii',Flag_frame);
    rv_TCR_c(ii,:) = [r_TCR_c_ii',v_TCR_c_ii'];
end

end

function [r_TCR_c,v_TCR_c] = TCO2TCR(r_TCO_c,v_TCO_c,r_MCR_t,v_MCR_t,a_MCR_t,Flag_frame)

h_MCR_t = cross(r_MCR_t,v_MCR_t);
h_MCR_t_norm = norm(h_MCR_t);
v_MCR_t_norm = norm(v_MCR_t);
r_MCR_t_norm = norm(r_MCR_t);
% ����ϵ����
if strcmp(Flag_frame,'VNC')
    iVector = v_MCR_t/v_MCR_t_norm;
    jVector = h_MCR_t/h_MCR_t_norm;
    kVector = cross(iVector,jVector);
    omegaT = (dot(a_MCR_t,h_MCR_t)*cross(v_MCR_t,h_MCR_t)...
        + dot(a_MCR_t,h_MCR_t)*dot(r_MCR_t,v_MCR_t)*v_MCR_t...
        - dot(a_MCR_t,cross(v_MCR_t,h_MCR_t))*h_MCR_t)/h_MCR_t_norm^2/v_MCR_t_norm^2;
elseif strcmp(Flag_frame,'LVLH')
    iVector = r_MCR_t/r_MCR_t_norm;
    kVector = h_MCR_t/h_MCR_t_norm;
    jVector = cross(kVector,iVector);
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
%     omegaT = h_MCR_t/r_MCR_t_norm^2;
elseif strcmp(Flag_frame,'VVLH')
    kVector = -r_MCR_t/r_MCR_t_norm;
    jVector = -h_MCR_t/h_MCR_t_norm;
    iVector = cross(jVector,kVector);
    omegaT = h_MCR_t/r_MCR_t_norm^2 + dot(a_MCR_t,h_MCR_t)/h_MCR_t_norm^2*r_MCR_t;
else
    error('Wrong orbital frame flag')
end

%����ת��
TM_TCO2TCR = [iVector jVector kVector];
r_TCR_c = TM_TCO2TCR*r_TCO_c;
% r_TCR_c = r_TCR_c;
% rMC_TCR_c = r_TCR_t+rTCO_TCR;

if isempty(v_TCO_c)
    v_TCR_c = [];
else
    v_TCR_c = TM_TCO2TCR*v_TCO_c + cross(omegaT,r_TCR_c);
end
% v_TCR_c = v_TCR_c;
% vMC_TCR_c = v_TCR_t+vTCO_TCR;

end

