function aa_j2k = eomj2kMtx(xx_j2k,t_sec,aux)
%% 计算j2k坐标系下的轨道加速度，以便相对运动坐标系的轨道转移使用
aa_j2k = zeros(size(xx_j2k));
if isfield(aux,'jd0_2et')
    for ii_loop = 1:length(t_sec)
        aa_j2k(ii_loop,:) = eom_SPICE_stm(t_sec(ii_loop), xx_j2k(ii_loop,:), aux);
    end
else
    for ii_loop = 1:length(t_sec)
        aa_j2k(ii_loop,:) = eqom_geoMEMEJ2k(t_sec(ii_loop), xx_j2k(ii_loop,:), aux);
    end
end
aa_j2k = aa_j2k(:,4:6);
