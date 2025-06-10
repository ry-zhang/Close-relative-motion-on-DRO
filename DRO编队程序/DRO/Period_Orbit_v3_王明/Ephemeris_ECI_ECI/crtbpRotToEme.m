
function rv_N3_eme = crtbpRotToEme(jdEpoch , P1 , P2 , centerflag , rv_O3_rot , aux)
%
% rot -> eme(星历)
%
% 输入:
% jdEpoch            [1x1]         历元
% P1                     [1x1]         P1编号
% P2                     [1x1]         P2编号
% centerflag         [1x1]         中心编号
% rv_O3_rot          [6x1]         rot状态
% aux                    [1x1]         参数
%
% 输出:
% rv_N3_eme             [6x1]         eme状态
%
% 参考：
% [1]: 杨洪伟，全星历模型下拟Halo轨道设计，2015
% [2]: 李明涛，共线平动点任务节能轨道设计与优化，2010
%
% 注意：
% ref [1]的omega计算错误
%
% 作者: 张晨, 中科院空间应用工程与技术中心
% chenzhang@csu.ac.cn
% 2019/12/26
% 2020/06/20
% --------------------------------------------------------------

rv_O3_rot = rv_O3_rot(:);

% 计算旋转系引力常数
mu = aux.planet.mu(P2) / (aux.planet.mu(P1) + aux.planet.mu(P2));
% mu = aux.mu;

% 旋转系参数
r_O1_rot = [-mu ; 0 ; 0];
r_O2_rot = [1 - mu ; 0 ; 0];
r_O3_rot = rv_O3_rot(1:3);
v_O3_rot = rv_O3_rot(4:6);

% -------------------------------------------------
% P2相对P1状态(eme)
rv_12_eme = ephEme(jdEpoch, P2 , P1 , aux.DE430);
r_12_eme = rv_12_eme(1:3);
v_12_eme = rv_12_eme(4:6);

% 计算瞬时角速率
omega_eme = cross(r_12_eme , v_12_eme) / norm(r_12_eme)^2;

% % 瞬时距离
LU = norm(r_12_eme);

% 瞬时速度
VU = norm(v_12_eme);

% 计算旋转矩阵(rot -> eme)
e1 = r_12_eme / norm(r_12_eme);
e3 = omega_eme / norm(omega_eme);
e2 = cross(e3 , e1);
M_rotToEme = [e1 , e2 , e3];

switch centerflag
    
    case 'P1'
        
        % 计算3相对1的位置矢量在eme投影
        r_23_eme = M_rotToEme * (r_O3_rot - r_O2_rot) * LU;
        r_13_eme = r_12_eme + r_23_eme;
        
        % 计算3相对1的速度矢量在eme投影
        v_13_eme = M_rotToEme * v_O3_rot * VU + cross(omega_eme , r_23_eme) + v_12_eme;
        rv_13_eme = [r_13_eme ; v_13_eme];
        
        rv_N3_eme = rv_13_eme;
        
    case 'P2'
        
        % 计算3相对2的位置矢量在eme投影
        r_13_eme = M_rotToEme * (r_O3_rot - r_O1_rot) * LU;
        r_23_eme = r_13_eme - r_12_eme;
        
        % 计算3相对1的速度矢量在eme投影
        v_23_eme = M_rotToEme * v_O3_rot * VU + cross(omega_eme , r_13_eme) - v_12_eme;
        rv_23_eme = [r_23_eme ; v_23_eme];
        
        rv_N3_eme = rv_23_eme;
        
end

end
