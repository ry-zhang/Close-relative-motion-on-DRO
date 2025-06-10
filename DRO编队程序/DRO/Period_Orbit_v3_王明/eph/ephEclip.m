
function rv = ephEclip(jd , ntarg, ncent , C_Mat)
%
% ¼ÆËãeclipÐÇÀú
%
%%%%%%%%%%%%%%%%%%%%%%%

% eme2000 -> eclip
eq2000 = [[ +1.000000000000d0 -0.000000479966d0 0.000000000000d0]; ...
    [ +0.000000440360d0 +0.917482137087d0 +0.397776982902d0]; ...
    [ -0.000000190919d0 -0.397776982902d0 +0.917482137087d0]];

rvCent = jplEph(jd, ncent, C_Mat);
rvTarg = jplEph(jd, ntarg, C_Mat);

rv = [eq2000 * (rvTarg(:, 1) - rvCent(: , 1));
    eq2000 * (rvTarg(:, 2) - rvCent(: , 2)) / 86400];

end
