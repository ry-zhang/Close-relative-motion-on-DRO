function [xx , tt] = crtbpP2centered2synodic3D(XX , TT , mu)
% Coordinate change from inertial, P2-centerd to synodic frame
%
% [xx, tt] = P2centered2synodic3D (XX, TT, mu)
%
% Input arguments:
% -------------------------------------------------------------------------
% XX        [nx6]    state vector in P2-centerd frame    [adim.]
% TT        [nx1]    time vector                         [adim.]
% mu        [1x1]    mass parameter of the RTBP          [adim.]
%
% Output arguments:
% -------------------------------------------------------------------------
% xx        [nx6]     state vector in synodic frame      [adim.]
% tt        [nx1]     time vector in synodic frame       [adim.]
%
% External functions called:
% -------------------------------------------------------------------------
% no calls
%
% Copyright (C) 2020/02/08 by Chen Zhang
sin_TT = sin(TT);
cos_TT = cos(TT);
xx = [XX(: , 1) .* cos_TT - mu + XX(: , 2) .* sin_TT + 1, ...
    XX(: , 2) .* cos_TT - XX(: , 1) .* sin_TT, ...
    XX(: , 3), ...
    XX(: , 4) .* cos_TT + XX(: , 2) .* cos_TT + XX(: , 5) .* sin_TT - XX(: , 1) .* sin_TT, ...
    XX(: , 5) .* cos_TT - XX(: , 1) .* cos_TT - XX(: , 4) .* sin_TT - XX(: , 2) .* sin_TT, ...
    XX(: , 6)];
tt = TT;
end