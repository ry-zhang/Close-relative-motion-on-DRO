function xx0 = FirstGuessFixedAz3D (Az, L, mu, PointFlag,Orbit_class)

% First guess state for Halo orbits, fixed Az
%
% xx0 = FirstGuessFixedAz3D (Az, L, mu, PointFlag)
%
% Input arguments:
% -------------------------------------------------------------------------
% Az         [1x1]    initial 3D amplitude                     [adim., EM units]
% L          [1x1]    initial in-plane amplitude               [adim.]
%                     L = L1   ->  L1 position                 [adim.]
%                     L = L2   ->  L2 position                 [adim.]
% mu         [1x1]    mass parameter of the P1-P2 r3bp         [adim.]
% PointFlag  [1x1]    1 -> L1, 2 -> L2                         [adim.]

%
% Output arguments:
% -------------------------------------------------------------------------
% xx0        [6x1]   first guess state xx0 = [x0 0 z0 0 vy0 0]'[adim.]
%
% External functions called:
% -------------------------------------------------------------------------
% none
%
% Copyright (C) 27/03/2012 by Renyong Zhang


% The coefficients for the definition of the first guess solution are taken by
% R. Thurman and P. Worfolk, The Geometry of Halo Orbits in the Circular Restricted
% Three-Body Problem, Geometry Center Research Report, University of Minnesota, GCG95, 1996.

% c2, c3, c4
c2=mu/(abs(L-1+mu))^3+(1-mu)/(abs(L+mu))^3 ;
if strcmp(PointFlag,'L1') % L1
    c3 = mu/(abs(L-1+mu))^3 - (1-mu)*abs(L-1+mu)/(abs(L+mu))^4 ;
elseif strcmp(PointFlag,'L2')
    c3 = -mu/(abs(L-1+mu))^3 - (1-mu)*abs(L-1+mu)/(abs(L+mu))^4 ;
end
c4 = mu/(abs(L-1+mu))^3 + (1-mu)*abs(L-1+mu)^2/(abs(L+mu))^5 ;

% lam, k, Delta
lam = sqrt((2-c2+sqrt(9*c2^2-8*c2))/2);
k = 2*lam/(lam^2+1-c2);
Delta = lam^2-c2;

% d1, d2, d3
d1 = 16*lam^4+4*lam^2*(c2-2)-2*c2^2+c2+1;
d2 = 81*lam^4+9*lam^2*(c2-2)-2*c2^2+c2+1;
d3 = 2*lam*(lam*(1+k^2)-2*k);

% aij, bij, dij
a21 = 3*c3*(k^2-2)/(4*(1+2*c2));
a22 = 3*c3/(4*(1+2*c2));
a23 = -3*lam*c3/(4*k*d1)*(3*k^3*lam-6*k*(k-lam)+4);
a24 = -3*lam*c3/(4*k*d1)*(2+3*lam*k);
b21 =-3*c3*lam/(2*d1)*(3*lam*k-4);
b22 = 3*lam*c3/d1;
d21 = -c3/(2*lam^2);
a31 = -9*lam/d2*(c3*(k*a23-b21) + k*c4*(1+1/4*k^2)) + (9*lam^2+1-c2)/(2*d2)*(3*c3*(2*a23-k*b21) + c4*(2+3*k^2));
a32 = -9*lam/(4*d2)*(4*c3*(k*a24-b22)+k*c4)-3*(9*lam^2+1-c2)/(2*d2)*(c3*(k*b22+d21-2*a24)-c4);
b31 = 1/d2*(3*lam*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2)) + (9*lam^2+1+2*c2)*(12*c3*(k*a23-b21) + 3*k*c4*(4+k^2))/8);
b32 = 1/d2*(3*lam*(3*c3*(k*b22+d21-2*a24)-3*c4) + (9*lam^2+1+2*c2)*(12*c3*(k*a24-b22) + 3*c4*k)/8);
d31 = 3/(64*lam^2)*(4*c3*a24+c4);
d32 = 3/(64*lam^2)*(4*c3*(a23-d21) + c4*(4+k^2));

% s1, s2
s1 = 1/d3*(3/2*c3*(2*a21*(k^2-2)-a23*(k^2+2)-2*k*b21)-3/8*c4*(3*k^4-8*k^2+8));
s2 = 1/d3*(3/2*c3*(2*a22*(k^2-2)+a24*(k^2+2)+2*k*b22+d21*(2+3))+3/8*c4*((8+4)-k^2*(2-1)));

% Parameters for the Thurman and Worfolk correction
b33 = -k/(16*lam)*(12*c3*(b21-2*k*a21+k*a23)+3*c4*k*(3*k^2-4)+16*s1*lam*(lam*k-1));
b34 = -k/(8*lam)*(-12*c3*k*a22+3*c4*k+8*s2*lam*(lam*k-1));
b35 = -k/(16*lam)*(12*c3*(b22+k*a24)+3*c4*k);

% a1, a2, l1, l2
a1 = -3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-k^2);
a2 = 3/2*c3*(a24-2*a22)+9/8*c4;
l1 = a1+2*lam^2*s1;
l2 = a2+2*lam^2*s2;

% scale Az (reference frame centered at L with length unit = L-P2 distance)
Az = Az/abs(L-1+mu) ;

% Az, Ax, Ay (adimensional)
Ax = sqrt(-(l2*Az^2+Delta)/l1); % Ax=Ax(Az) subject to a constraint !!!
Ay = k*Ax;

% First-guess initial condition (third order Lindstedt-Poincar?method)
% x0
x0 = -Ax + a21*Ax^2 + a22*Az^2 + a23*Ax^2 - a24*Az^2 + a31*Ax^3 - a32*Ax*Az^2;

% z0
z0 = Az - 2*d21*Ax*Az + d32*Az*Ax^2 - d31*Az^3;

% v0
om2 = s1*Ax^2 + s2*Az^2;
dtao1 = lam*(1+om2);
v0 = k*Ax*dtao1 + 2*dtao1*(b21*Ax^2-b22*Az^2) + 3*dtao1*(b31*Ax^3 - b32*Ax*Az^2);

% Thurman and Worfolk correction
v0_mod = v0 + dtao1*(b33*Ax^3 + b34*Ax*Az^2 - b35*Ax*Az^2); 

% re-transform x0 and vy0 in the RTBP synodic frame
x0 = L + x0*abs(L-1+mu) ;
z0 = z0*abs(L-1+mu) ;
v0_mod = v0_mod*abs(L-1+mu) ;

%first-guessed xx0
if strcmp(PointFlag,'L1')
    if strcmp(Orbit_class,'North')
        xx0 = [x0 0 z0 0 v0_mod 0]' ;
    elseif strcmp(Orbit_class,'South')
        xx0 = [x0 0 -z0 0 v0_mod 0]' ;
    end
    
elseif strcmp(PointFlag,'L2')
    if strcmp(Orbit_class,'North')
        xx0 = [x0 0 -z0 0 v0_mod 0]' ;
    elseif strcmp(Orbit_class,'South')
        xx0 = [x0 0 z0 0 v0_mod 0]' ;
    end
end
end

