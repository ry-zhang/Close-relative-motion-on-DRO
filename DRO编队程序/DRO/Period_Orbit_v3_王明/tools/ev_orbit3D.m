function [value,isterminal,direction] = ev_orbit3D (t, XM, aux)

% Event Function for Initial Condition Orbit Correction
%
% [value,isterminal,direction] = ev_orbit (t, XM, mu)
%
% Input arguments:
% -------------------------------------------------------------------------
% t          [1x1]    time                                         [adim.]
% XM         [1x42]   6-state vector (adim) + 6x6 monodromy matrix [adim.]
% mu         [1x1]    mass parameter (adim)                        [adim.]
%
% Output arguments:
% -------------------------------------------------------------------------
% value      [1x1]     y-coordinate                                [adim.]
% isterminal [1x1]     stops at the first zero of "value" (1)      [adim.]
% direction  [1x1]     stops when "value" decrease (-1)            [adim.]
%
% External functions called:
% -------------------------------------------------------------------------
% none
%
% Copyright (C) 26/03/2012 by Renyong Zhang


value=XM(2);

isterminal=1;
direction=-1;

end