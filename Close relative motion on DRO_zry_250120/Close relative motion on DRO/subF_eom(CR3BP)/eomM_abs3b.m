function ydot = eomM_abs3b(t,y,mu)
% equation of motion of CR3BP, with equation of state of transition matrix
%   in synodic coordiante system cetered in Moon
%
% ydot = eomM_abs3b(t,y,mu)
%
% Input arguments:
% -------------------------------------------------------------------------
% t          [1x1]               FFT size
% y          [(6+36)x1]          
%                1:3             Position
%                4:6             Velocity
%                7:42            State transition matrix
% mu         [1x1]               
% 
% Output arguments:
% -------------------------------------------------------------------------
% ydot       [(6+36)x1]          dot Matrix
% 
% External functions called:
% -------------------------------------------------------------------------
%  none
% 
% Copyright (C) 25/7/2020 by Chihang Yang 
% email: ychhtl@foxmail.com
% v1(2022-03-01):
% ����ת����ϵ�ɵ������ұ��޸�Ϊ���������
%% ����/��ƹ��������������ϵM�µĶ���ѧ
% ԭ��λ��������ת����ϵM��x��ӵ���ָ������y�������Ƶ�����ת�ķ���
r_C = y(1:3); % ����/��ƹ����������ת����ϵM�µ�λ���ٶ�
v_C = y(4:6);
r_C_dot = v_C;

% ��ǰ��״̬ת�ƾ���
M = reshape(y(7:42),6,6);

% Բ������������
% �������������ȡ���ڻ������ϵ��ԭ���Ƿ������Ļ��Ƿ��������
r_E_C = r_C-[-1;0;0]; r_E_C_norm  = norm(r_E_C); 
r_E_C3 = r_E_C_norm^3; r_E_C5 = r_E_C_norm^5;
r_M_C_norm  = norm(r_C); r_M_C3 = r_M_C_norm^3; r_M_C5 = r_M_C_norm^5; 

% ����ѧ����
Ar_C = [1-mu/r_M_C3-(1-mu)/r_E_C3, 0, 0;
        0, 1-mu/r_M_C3-(1-mu)/r_E_C3, 0;
        0, 0, -mu/r_M_C3-(1-mu)/r_E_C3];
Av_C = [0, 2, 0;
        -2, 0, 0;
        0, 0, 0];
B = [(1-mu)*(1-1/r_E_C3); 0; 0];
v_C_dot = Ar_C*r_C + Av_C*v_C + B;
rv_dot = [r_C_dot; v_C_dot];

% ��״̬ת�ƾ���ĵ��������Խ��ƣ�
rx = r_C(1); ry = r_C(2); rz = r_C(3);
omegaxx = 1 + 3*mu*rx^2/r_M_C5 - mu/r_M_C3 - (1-mu)*(1/r_E_C3 - 3*(rx+1)^2/r_E_C5);
omegaxy = 3*(1-mu)*(rx+1)*ry/r_E_C5 + 3*mu*rx*ry/r_M_C5 ;
omegaxz = 3*(1-mu)*(rx+1)*rz/r_E_C5 + 3*mu*rx*rz/r_M_C5 ;
omegayy = 1 - (1-mu)/r_E_C3 + 3*(1-mu)*ry^2/r_E_C5 - mu/r_M_C3 + 3*mu*ry^2/r_M_C5 ;
omegayz = 3*(1-mu)*ry*rz/r_E_C5 + 3*mu*ry*rz/r_M_C5 ;
omegazz = -(1-mu)/r_E_C3 + 3*(1-mu)*rz^2/r_E_C5 - mu/r_M_C3 + 3*mu*rz^2/r_M_C5 ;


A = [0       0       0        1   0  0 ;
     0       0       0        0   1  0 ;
     0       0       0        0   0  1 ;
     omegaxx omegaxy omegaxz  0   2  0 ;
     omegaxy omegayy omegayz  -2  0  0 ;
     omegaxz omegayz omegazz  0   0  0] ;

Mdot = A*M ;

ydot = [rv_dot; Mdot(:)];
