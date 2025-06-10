function [rv_DRO_MCR,a_MCR] = ECJ2k2Rot(t_sample,rv_DRO,aux)
size_Data = length(t_sample);
a_all = zeros(size_Data,3);
parfor ii_loop = 1:size_Data
    x_temp = rv_DRO(ii_loop,:);
    t_temp = t_sample(ii_loop);
    a_temp = eqom_geoMEMEJ2k(t_temp, x_temp, aux);
    a_all(ii_loop,:) = a_temp(4:6);
end
tt_jd = aux.jd0 + t_sample/86400;
[rv_DRO_MCR,a_MCR] = T_ECJ2k2Rot(tt_jd, rv_DRO, a_all,aux.C_Mat, 'MCEMR');
