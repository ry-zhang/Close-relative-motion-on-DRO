function x0 = Refine_x0(x0, TT, aux)


Tspan = [0, TT];
options = odeset('RelTol',1e-13,'AbsTol',1e-16);
[tt, xx]  = ode113(@crtbpEqm3D, Tspan, x0 , options, aux);
% 三次样条插值
if strcmp(aux.periOrb_Modi(1:2) , 'RO')
    num_str = regexp(aux.periOrb_Modi,'\d*\.?\d*','match');
    num = str2double(num_str);
    
    periOrb_p = num(1);
    aux.periOrb_n = aux.periOrb_n*periOrb_p;
end
tt_interp = linspace(0 , TT , aux.periOrb_n)';
xx_interp = interp1(tt , xx , tt_interp' , 'spline')' ;
xj = xx_interp(:);
err = inf;
k = 0;
aux.yFix = x0(2);
aux.Period = TT;
while err > 1e-10  
    aux.pseudoArc_IO = 0;
    
    [~, cc_eq, ~, Gc_eq] = constr_cr3bp(xj, aux);
    Gc_eq = Gc_eq';
    
    % 计算误差
    err = max(abs(cc_eq));
    fprintf('err: %0.4e \n' , err)
    
    % 更新xj
    xj = xj - Gc_eq' * pinv(Gc_eq * Gc_eq') * cc_eq; 
    k = k+1;
end
    x0 = xj(1:6);
end

