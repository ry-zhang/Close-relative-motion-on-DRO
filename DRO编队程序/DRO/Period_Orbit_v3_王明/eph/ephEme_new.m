
function rv = ephEme_new(jd , ntarg, ncent , C_Mat)
%
% ¼ÆËãeme2000ÐÇÀú
% 
%%%%%%%%%%%%%%%%%%%%%%%

rvCent = jplEph_new(jd, ncent, C_Mat);
rvTarg = jplEph_new(jd, ntarg, C_Mat);

rv = [rvTarg(:, 1) - rvCent(: , 1);
    (rvTarg(:, 2) - rvCent(: , 2)) / 86400;
    (rvTarg(:, 3) - rvCent(: , 3)) / 86400^2;
    (rvTarg(:, 4) - rvCent(: , 4)) / 86400^3];

end
