function [position,isterminal,direction] = DRO_eph_event2(t, y, aux, target_theta)
x_j2k = y';

% DRO ECIj2k → MCR
tt_jd = aux.jd0 + t/86400;
x_MCR = T_ECJ2k2Rot(tt_jd, x_j2k, [],aux.C_Mat, 'MCEMR',0);

position = (mod(atan2(-x_MCR(2),x_MCR(1)),2*pi) - target_theta); % The value that we want to be zero
isterminal = zeros(size(target_theta(:)));  % Halt integration 
direction = ones(size(target_theta(:)));   % Positive direction only