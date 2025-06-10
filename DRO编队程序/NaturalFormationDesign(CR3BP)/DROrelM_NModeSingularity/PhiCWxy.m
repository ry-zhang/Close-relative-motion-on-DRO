function Phi_t = PhiCWxy(t,n)
nt = n*t;
ct = cos(nt);
st = sin(nt);
Phi_t = [4-3*ct 0 st/n 2/n-2*ct/n;
    -6*nt+6*st 1 -2/n+2*ct/n 4*st/n-3*t;
    3*n*st 0 ct 2*st;
    -6*n+6*n*ct 0 -2*st -3+4*ct];