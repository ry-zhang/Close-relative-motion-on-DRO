function [L1, L2, L3, L4x, L4y] = Li (mu)

% Lagrangian points location of the RTBP
%
% [L1, L2, L3, L4x, L4y] = Li (mu)
%
% Input argumuents:
% -------------------------------------------------------------------------
% mu        [1x1]     mass parameter of the r3bp                  [adim.]
%
% Output argumuents:
% -------------------------------------------------------------------------
% L1        [1x1]     x-ccordinate of L1                          [adim.]
% L2        [1x1]     x-ccordinate of L2                          [adim.]
% L3        [1x1]     x-ccordinate of L3                          [adim.]
% L4x       [1x1]     x-coordinate of L4 (=L5x)                   [adim.]
% L4y       [1x1]     y-coordinate of L5 (=-L5y)                  [adim.]
% -------------------------------------------------------------------------
%
% External functions called:
% -------------------------------------------------------------------------
% none
%
% Copyright (C) 15/7/2012 by by Renyong Zhang
% Modified on 07/03/2012 by Renyong Zhang

options = optimset('TolFun', 1e-14, 'TolX', 1e-14) ;

x012 = (mu/3)^(1/3) ;
x03  = 1-(7/12)*mu  ;

L1  = 1 - mu - fzero(@L1poly, x012, options, mu) ;
L2  = 1 - mu + fzero(@L2poly, x012, options, mu) ;
L3  =   - mu - fzero(@L3poly, x03,  options, mu) ;
L4x = 0.5 - mu ;
L4y = sqrt(3)/2 ;

end

%

function [poly] = L1poly (x,mu)

poly=x.^5-(3-mu)*x.^4+(3-2*mu)*x.^3-mu*x.^2+2*mu*x-mu;

end

%

function [poly] = L2poly (x,mu)

poly=x.^5+(3-mu)*x.^4+(3-2*mu)*x.^3-mu*x.^2-2*mu*x-mu;

end

%

function [poly] = L3poly (x,mu)

poly=x.^5+(2+mu)*x.^4+(1+2*mu)*x.^3-(1-mu)*x.^2-2*(1-mu)*x-(1-mu);

end