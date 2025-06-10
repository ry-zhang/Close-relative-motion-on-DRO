
function rv_N3_eme = crtbpRotToEme(jdEpoch , P1 , P2 , centerflag , rv_O3_rot , aux)
%
% rot -> eme(����)
%
% ����:
% jdEpoch            [1x1]         ��Ԫ
% P1                     [1x1]         P1���
% P2                     [1x1]         P2���
% centerflag         [1x1]         ���ı��
% rv_O3_rot          [6x1]         rot״̬
% aux                    [1x1]         ����
%
% ���:
% rv_N3_eme             [6x1]         eme״̬
%
% �ο���
% [1]: ���ΰ��ȫ����ģ������Halo�����ƣ�2015
% [2]: �����Σ�����ƽ����������ܹ��������Ż���2010
%
% ע�⣺
% ref [1]��omega�������
%
% ����: �ų�, �п�Ժ�ռ�Ӧ�ù����뼼������
% chenzhang@csu.ac.cn
% 2019/12/26
% 2020/06/20
% --------------------------------------------------------------

rv_O3_rot = rv_O3_rot(:);

% ������תϵ��������
mu = aux.planet.mu(P2) / (aux.planet.mu(P1) + aux.planet.mu(P2));
% mu = aux.mu;

% ��תϵ����
r_O1_rot = [-mu ; 0 ; 0];
r_O2_rot = [1 - mu ; 0 ; 0];
r_O3_rot = rv_O3_rot(1:3);
v_O3_rot = rv_O3_rot(4:6);

% -------------------------------------------------
% P2���P1״̬(eme)
rv_12_eme = ephEme(jdEpoch, P2 , P1 , aux.DE430);
r_12_eme = rv_12_eme(1:3);
v_12_eme = rv_12_eme(4:6);

% ����˲ʱ������
omega_eme = cross(r_12_eme , v_12_eme) / norm(r_12_eme)^2;

% % ˲ʱ����
LU = norm(r_12_eme);

% ˲ʱ�ٶ�
VU = norm(v_12_eme);

% ������ת����(rot -> eme)
e1 = r_12_eme / norm(r_12_eme);
e3 = omega_eme / norm(omega_eme);
e2 = cross(e3 , e1);
M_rotToEme = [e1 , e2 , e3];

switch centerflag
    
    case 'P1'
        
        % ����3���1��λ��ʸ����emeͶӰ
        r_23_eme = M_rotToEme * (r_O3_rot - r_O2_rot) * LU;
        r_13_eme = r_12_eme + r_23_eme;
        
        % ����3���1���ٶ�ʸ����emeͶӰ
        v_13_eme = M_rotToEme * v_O3_rot * VU + cross(omega_eme , r_23_eme) + v_12_eme;
        rv_13_eme = [r_13_eme ; v_13_eme];
        
        rv_N3_eme = rv_13_eme;
        
    case 'P2'
        
        % ����3���2��λ��ʸ����emeͶӰ
        r_13_eme = M_rotToEme * (r_O3_rot - r_O1_rot) * LU;
        r_23_eme = r_13_eme - r_12_eme;
        
        % ����3���1���ٶ�ʸ����emeͶӰ
        v_23_eme = M_rotToEme * v_O3_rot * VU + cross(omega_eme , r_13_eme) - v_12_eme;
        rv_23_eme = [r_23_eme ; v_23_eme];
        
        rv_N3_eme = rv_23_eme;
        
end

end
