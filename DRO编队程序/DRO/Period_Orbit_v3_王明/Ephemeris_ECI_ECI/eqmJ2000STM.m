
function dxx = eqmJ2000STM(t , xx , aux)
%
% J2000动力学方程
%
%               i (SC)
%             /    \
%           /          \
%         /               \
%   k(E) -------------- j (M/S)
%
% 参考：
% [1]: Vaquero. SpaceCraft Transfer Trajectory Design Exploiting Resonant Orbits in
% Multi-Body Environments [D], Purdue University, 2013, page 33-35.
% [2]: Thomas A. Pavlak. Trajectory design and Orbit Maintenance Strategies in Multi-Body
% Dynamical Regimes [D], Purdue University, 2013, page 34-41.
% 
% 输入：
% t: jdtdb(s)
% xx: j2000状态(km, km/s)
% k_index: 主天体编号
% j_index_vec: 摄动天体编号
%
% 作者: 张晨, 中科院空间应用工程与技术中心
% chenzhang@csu.ac.cn
% 2021/06/06
% -----------------------------------------------------------

% 行向量 -> 列向量
xx = xx(:);

% 出发时间
dep_t_jdtdb = aux.dep_t_jdtdb;

k_index = 3;
j_index = [10 , 11];

if size(xx , 1) == 6
    STM_IO = 0;
else
    STM_IO = 1;
end

% ------------------------ k planet ------------------------
% r_ki / R_ki
r_ki = xx(1 : 3);
R_ki = sqrt(r_ki(1)^2 + r_ki(2)^2 + r_ki(3)^2);
R_ki_3 = R_ki^3;
R_ki_5 = R_ki^5;

% update vector_field
mu_k = aux.planet.mu(k_index) * aux.TU^2 / aux.LU^3;
acc = - mu_k * r_ki / R_ki_3;

% update STM
if STM_IO == 1
    
    A41 = - mu_k * (1 / R_ki_3 - 3 * r_ki(1) * r_ki(1) / R_ki_5);
    A52 = - mu_k * (1 / R_ki_3 - 3 * r_ki(2) * r_ki(2) / R_ki_5);
    A63 = - mu_k * (1 / R_ki_3 - 3 * r_ki(3) * r_ki(3) / R_ki_5);
    
    A42 = - mu_k * (-3 * r_ki(1) * r_ki(2) / R_ki_5);
    A43 = - mu_k * (-3 * r_ki(1) * r_ki(3) / R_ki_5);
    A53 = - mu_k * (-3 * r_ki(2) * r_ki(3) / R_ki_5);
    
    A51 = A42;
    A61 = A43;
    A62 = A53;
    
end

TEMP = zeros(6 , 1);
% ------------------------ j planet ------------------------
for iLoop = 1 : length(j_index)
    
    % r_kj
    x_temp = ephEme_mex(dep_t_jdtdb + t * aux.TU / 86400 , j_index(iLoop) , k_index , aux.DE430);
    rv_kj = [x_temp(1:3) / aux.LU; x_temp(4:6) / aux.VU];
    
    r_kj = rv_kj(1:3);
    v_kj = rv_kj(4:6);
    
    R_kj = sqrt(r_kj(1)^2 + r_kj(2)^2 + r_kj(3)^2);
    R_kj_3 = R_kj^3;
    R_kj_5 = R_kj^5;
    
    % r_ji
    r_ji = r_ki - r_kj;
    R_ji = sqrt(r_ji(1)^2 + r_ji(2)^2 + r_ji(3)^2);
    R_ji_3 = R_ji^3;
    R_ji_5 = R_ji^5;
    
    % acc_j
    mu_j = aux.planet.mu(j_index(iLoop)) * aux.TU^2 / aux.LU^3;
    acc = acc - mu_j * (r_kj / R_kj_3 + r_ji / R_ji_3);
%     acc = acc/aux.LU^2;
    B = zeros(6,3);
    B(4,1) = - mu_j * (1 / R_kj_3 - 3 * r_kj(1) * r_kj(1) / R_kj_5 - 1 / R_ji_3 + 3 * (r_ki(1) - r_kj(1)) * (r_ki(1) - r_kj(1)) / R_ji_5);
    B(5,2) = - mu_j * (1 / R_kj_3 - 3 * r_kj(2) * r_kj(2) / R_kj_5 - 1 / R_ji_3 + 3 * (r_ki(2) - r_kj(2)) * (r_ki(2) - r_kj(2)) / R_ji_5);
    B(6,3) = - mu_j * (1 / R_kj_3 - 3 * r_kj(3) * r_kj(3) / R_kj_5 - 1 / R_ji_3 + 3 * (r_ki(3) - r_kj(3)) * (r_ki(3) - r_kj(3)) / R_ji_5);
    B(4,2) = - mu_j * (- 3 * r_kj(1) * r_kj(2) / R_kj_5 + 3 * (r_ki(1) - r_kj(1)) * (r_ki(2) - r_kj(2)) / R_ji_5);
    B(4,3) = - mu_j * (- 3 * r_kj(1) * r_kj(3) / R_kj_5 + 3 * (r_ki(1) - r_kj(1)) * (r_ki(3) - r_kj(3)) / R_ji_5);
    B(5,3) = - mu_j * (- 3 * r_kj(2) * r_kj(3) / R_kj_5 + 3 * (r_ki(2) - r_kj(2)) * (r_ki(3) - r_kj(3)) / R_ji_5);
    B(5,1) = B(4,2);
    B(6,1) = B(4,3);
    B(6,2) = B(5,3);
    TEMP = TEMP + B * v_kj;
    
    % update STM
    if STM_IO == 1
        
        A41 = A41 - mu_j * (1 / R_ji_3 - 3 * r_ji(1) * r_ji(1) / R_ji_5);
        A52 = A52 - mu_j * (1 / R_ji_3 - 3 * r_ji(2) * r_ji(2) / R_ji_5);
        A63 = A63 - mu_j * (1 / R_ji_3 - 3 * r_ji(3) * r_ji(3) / R_ji_5);
        
        A42 = A42 - mu_j * (-3 * r_ji(1) * r_ji(2) / R_ji_5);
        A43 = A43 - mu_j * (-3 * r_ji(1) * r_ji(3) / R_ji_5);
        A53 = A53 - mu_j * (-3 * r_ji(2) * r_ji(3) / R_ji_5);
        
        A51 = A42;
        A61 = A43;
        A62 = A53;
        
    end
    
end

% ------------------------- vector field of state & STM --------------------------
if STM_IO == 0
    
    dxx = [xx(4:6); acc];
    
else
    
    % STM
    Phi = reshape(xx(7 : 42) , 6 , 6);
    
    % vector field of STM
    AMtx = [0 , 0 , 0 , 1 , 0 , 0;
        0 , 0 , 0 , 0 , 1 , 0;
        0 , 0 , 0 , 0 , 0 , 1;
        A41 , A42 , A43 , 0 , 0 , 0;
        A51 , A52 , A53 , 0 , 0 , 0;
        A61 , A62 , A63 , 0 , 0 , 0];
    dPhi = AMtx * Phi;
    
    dxx = [xx(4:6); acc; dPhi(:) ; AMtx * xx(43:48) + TEMP];
    
end

end
