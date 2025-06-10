
function rv_O3_rot = crtbpEmeToRot(jdEpoch , P1 , P2 , centerflag , rv_N3_eme , aux)
%
% eme -> rot(星历)
%
% 输入:
% jdEpoch            [1x1]         历元
% P1                     [1x1]         P1编号
% P2                     [1x1]         P2编号
% centerflag         [1x1]         中心编号
% PtoS_eme            [6x1]         eme状态
% aux                   [1x1]         参数
%
% 输出:
% OtoS_rot           [6x1]         rot状态
%
%
% 参考：
% [1]: 杨洪伟，全星历模型下拟Halo轨道设计，2015
% [2]: 李明涛，共线平动点任务节能轨道设计与优化，2010
%
% 作者: 张晨, 中科院空间应用工程与技术中心
% chenzhang@csu.ac.cn
% 2019/12/26
% --------------------------------------------------------------

rv_N3_eme = rv_N3_eme(:);

% 计算旋转系引力常数
mu = aux.planet.mu(P2) / (aux.planet.mu(P1) + aux.planet.mu(P2));

% 旋转系参数
r_O1_rot = [-mu ; 0 ; 0];
r_O2_rot = [1 - mu ; 0 ; 0];

% -------------------------------------------------
% P2相对P1状态(eme)
rv_12_eme = ephEme(jdEpoch , P2 , P1 , aux.DE430);
r_12_eme = rv_12_eme(1:3);
v_12_eme = rv_12_eme(4:6);

% 计算瞬时角速率
omega = cross(r_12_eme , v_12_eme) / norm(r_12_eme)^2;

% 瞬时距离
LU = norm(r_12_eme);

% 瞬时速度
VU = norm(omega) * LU;

% 计算旋转矩阵(rot -> eme)
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
        
        % 计算3相对O的位置矢量在rot投影
        r_O3_rot = M_emeToRot * (r_13_eme - r_12_eme) / LU + r_O2_rot;
        
        % 计算3相对O的速度矢量在rot投影
        v_O3_rot = M_emeToRot *  ( v_13_eme - v_12_eme - cross(omega , r_23_eme) ) / VU;
        
    case 'P2'
        
        r_23_eme = rv_N3_eme(1:3);
        v_23_eme = rv_N3_eme(4:6);
        r_13_eme = r_23_eme + r_12_eme;
        
        % 计算3相对O的位置矢量在rot投影
        r_O3_rot = M_emeToRot * (r_23_eme + r_12_eme) / LU + r_O1_rot;
        
        % 计算3相对O的速度矢量在rot投影
        v_O3_rot = M_emeToRot * ( v_23_eme + v_12_eme - cross(omega , r_13_eme) ) / VU;
        
end

rv_O3_rot = [r_O3_rot ; v_O3_rot];

end
