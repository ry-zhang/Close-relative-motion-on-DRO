function coeff = getCoeff(icStart, nc, ns, jIntv,PCtemp)
PCtemp = PCtemp(icStart : icStart - 1 + 3 * nc * ns); % all subintervals
PCtemp = reshape(PCtemp',  nc*3, ns);
PCtemp = PCtemp';
coeffTemp = PCtemp(jIntv , :); % get the correct subinterval
coeff = reshape(coeffTemp,nc,3);
coeff = coeff';
end