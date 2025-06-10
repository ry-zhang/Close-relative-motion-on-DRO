dv = zeros(5,4);
dv(1,1) = auxSFF.sepSFF.t0;
dv(2,1) = auxSFF.transferSFF(1).t0;
dv(3,1) = auxSFF.transferSFF(1).tf;
dv(4,1) = auxSFF.transferSFF(2).t0;
dv(5,1) = auxSFF.transferSFF(2).tf;

dv(1,2:4) = auxSFF.sepSFF.dv0_j2kLVLH;

dv(2,2:4) = auxSFF.transferSFF(1).dv0_j2kLVLH;
dv(3,2:4) = auxSFF.transferSFF(1).dvf_j2kLVLH;
dv(4,2:4) = auxSFF.transferSFF(2).dv0_j2kLVLH;
dv(5,2:4) = auxSFF.transferSFF(2).dvf_j2kLVLH;




save dv dv