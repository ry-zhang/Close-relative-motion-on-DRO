

function [tt , xx] = crtbpFlow3D(xx0 , tf , linetype , linewidth , aux)
%
% crtbp三维数值积分
%
% 例子：
% crtbpMarkEM
% aux = crtbpAuxEM;
% xx0 = [0.75 , 0 , 0 , 0.15 , 0 , 0];
% crtbpFlow3D(xx0 , 10 , aux.mu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('Reltol' , aux.tol , 'AbsTol' , aux.tol);
[tt , xx] = ode113(@fx , [0 , tf] , xx0 , options , aux);
plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , linetype , 'linewidth' , linewidth);

% 动力学方程
    function dxxdt = fx (t, xx, aux)
        mu = aux.mu;
        r1cube = ((xx(1) + mu)^2 + xx(2)^2 + xx(3)^2)^(1.5);
        r2cube = ((xx(1) - 1 + mu)^2 + xx(2)^2 + xx(3)^2)^(1.5);
        dxxdt = [xx(4);
            xx(5);
            xx(6);
            xx(1) + 2*xx(5) - (1 - mu) * (xx(1) + mu) / r1cube - mu*(xx(1) - 1 + mu) / r2cube;
            xx(2) - 2*xx(4) - (1 - mu) * xx(2) / r1cube - mu * xx(2) / r2cube;
            -(1 - mu) * xx(3) / r1cube - mu * xx(3) / r2cube];
    end

end

