
function rv_O3_rot = crtbpEmeToRot(jdEpoch , P1 , P2 , centerflag , rv_N3_eme , aux)
%
% eme -> rot(����)
%
% ����:
% jdEpoch            [1x1]         ��Ԫ
% P1                     [1x1]         P1���
% P2                     [1x1]         P2���
% centerflag         [1x1]         ���ı��
% PtoS_eme            [6x1]         eme״̬
% aux                   [1x1]         ����
%
% ���:
% OtoS_rot           [6x1]         rot״̬
%
%
% �ο���
% [1]: ���ΰ��ȫ����ģ������Halo�����ƣ�2015
% [2]: �����Σ�����ƽ����������ܹ��������Ż���2010
%
% ����: �ų�, �п�Ժ�ռ�Ӧ�ù����뼼������
% chenzhang@csu.ac.cn
% 2019/12/26
% --------------------------------------------------------------

rv_N3_eme = rv_N3_eme(:);

% ������תϵ��������
mu = aux.planet.mu(P2) / (aux.planet.mu(P1) + aux.planet.mu(P2));

% ��תϵ����
r_O1_rot = [-mu ; 0 ; 0];
r_O2_rot = [1 - mu ; 0 ; 0];

% -------------------------------------------------
% P2���P1״̬(eme)
rv_12_eme = ephEme(jdEpoch , P2 , P1 , aux.DE430);
r_12_eme = rv_12_eme(1:3);
v_12_eme = rv_12_eme(4:6);

% ����˲ʱ������
omega = cross(r_12_eme , v_12_eme) / norm(r_12_eme)^2;

% ˲ʱ����
LU = norm(r_12_eme);

% ˲ʱ�ٶ�
VU = norm(omega) * LU;

% ������ת����(rot -> eme)
e1 = r_12_eme / norm(r_12_eme);
e3 = omega / norm(omega);
e2 = cross(e3 , e1);
M_rotToEme = [e1 , e2 , e3];
M_emeToRot = M_rotToEme';

switch centerflag
    
    case 'P1'
        
        r_13_eme = rv_N3_eme(1:3);
        v_13_eme = rv_N3_eme(4:6);
        r_23_eme = r_13_eme - r_12_eme;
        
        % ����3���O��λ��ʸ����rotͶӰ
        r_O3_rot = M_emeToRot * (r_13_eme - r_12_eme) / LU + r_O2_rot;
        
        % ����3���O���ٶ�ʸ����rotͶӰ
        v_O3_rot = M_emeToRot *  ( v_13_eme - v_12_eme - cross(omega , r_23_eme) ) / VU;
        
    case 'P2'
        
        r_23_eme = rv_N3_eme(1:3);
        v_23_eme = rv_N3_eme(4:6);
        r_13_eme = r_23_eme + r_12_eme;
        
        % ����3���O��λ��ʸ����rotͶӰ
        r_O3_rot = M_emeToRot * (r_23_eme + r_12_eme) / LU + r_O1_rot;
        
        % ����3���O���ٶ�ʸ����rotͶӰ
        v_O3_rot = M_emeToRot * ( v_23_eme + v_12_eme - cross(omega , r_13_eme) ) / VU;
        
end

rv_O3_rot = [r_O3_rot ; v_O3_rot];

end
