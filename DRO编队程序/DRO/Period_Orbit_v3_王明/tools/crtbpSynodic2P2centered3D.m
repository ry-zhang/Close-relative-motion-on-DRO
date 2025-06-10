function [XX , TT] = crtbpSynodic2P2centered3D(xx , tt , mu)
% Coordinate change from synodic to inertial, P2-centered frame
%
% [XX, TT] = synodic2P2centered3D (xx, tt, mu)
%
% Input arguments:
% -------------------------------------------------------------------------
% xx        [nx6]    RTBP state vector in synodic frame  [adim.]
% tt        [nx1]    time vector in synodic frame        [adim.]
% mu        [1x1]    mass parameter of the RTBP          [adim.]
%
% Output arguments:
% -------------------------------------------------------------------------
% XX        [nx6]    state vector in Earth-centerd frame [adim.]
% TT        [nx1]    time vector                         [adim.]
%
% External functions called:
% -------------------------------------------------------------------------
% none
%
% Copyright (C) 2020/02/08 by Chen Zhang
sin_tt = sin(tt);
cos_tt = cos(tt);
XX = [cos_tt .* (mu + xx(: , 1) - 1) - xx(: , 2) .* sin_tt, ...
    sin_tt .* (mu + xx(: , 1) - 1) + xx(: , 2) .* cos_tt, ...
    xx(: , 3), ...
    cos_tt .* (xx(: , 4) - xx(: , 2)) - sin_tt .* (mu + xx(: , 5) + xx(: , 1) - 1), ...
    cos_tt .* (mu + xx(: , 5) + xx(: , 1) - 1) + sin_tt .* (xx(: , 4) - xx(: , 2)), ...
    xx(: , 6)];
TT = tt;
end