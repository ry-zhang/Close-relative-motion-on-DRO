function rv = cheb_Interp_all(coeff, tau,nc)
% km, 32/ns-day
rv = zeros(3,4);
pc = zeros(1 , nc);
vc = pc;
pc(1) = 1;
pc(2) = tau;
vc(1) = 0;
vc(2) = 1;
ac(1) = 0;
ac(2) = 0;
ac(3) = 4;
ajc(1) = 0;
ajc(2) = 0;
ajc(3) = 0;
ajc(4) = 24;
for ii = 3 : nc
    pc(ii) = 2*tau*pc(ii-1)-pc(ii-2);
    vc(ii) = 2*pc(ii-1)+2*tau*vc(ii-1)-vc(ii-2);
end
for ii =4:nc
    ac(ii) = 4*vc(ii-1)+2*tau*ac(ii-1)-ac(ii-2);
end
for ii =5 :nc
    ajc(ii) = 6*ac(ii-1) + 2*tau*ajc(ii-1) - ajc(ii-2);
end
pc(3 : nc) = 2 * tau * pc(2 : (nc -1)) - pc(1 : (nc - 2));

rv(: , 1) = coeff * pc';
rv(: , 2) = coeff * vc';
rv(: , 3) = coeff * ac';
rv(: , 4) = coeff * ajc';
end
