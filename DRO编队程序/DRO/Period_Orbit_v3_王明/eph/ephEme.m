
function rv = ephEme(jd , ntarg, ncent , C_Mat)
%
% ����eme2000����
% 
%%%%%%%%%%%%%%%%%%%%%%%

rvCent = jplEph(jd, ncent, C_Mat);
rvTarg = jplEph(jd, ntarg, C_Mat);

rv = [rvTarg(:, 1) - rvCent(: , 1);
    (rvTarg(:, 2) - rvCent(: , 2)) / 86400];

end
