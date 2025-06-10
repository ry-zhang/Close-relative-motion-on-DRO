function rv = cheb_Interp(coeff, tau,nc)
% km, 32/ns-day
rv = zeros(3,2);
pc = zeros(1 , nc);
vc = pc;
pc(1) = 1;
pc(2) = tau;
vc(1) = 0;
vc(2) = 1;

for ii = 3 : nc
    pc(ii) = 2*tau*pc(ii-1)-pc(ii-2);
    vc(ii) = 2*pc(ii-1)+2*tau*vc(ii-1)-vc(ii-2);
end

pc(3 : nc) = 2 * tau * pc(2 : (nc -1)) - pc(1 : (nc - 2));

rv(: , 1) = coeff * pc';
rv(: , 2) = coeff * vc';

end
