function dydt = odeR3bpAug(t,  y , auxdata)
% 限制性三体问题
% 如果dim(y)>6, 计算状态转移矩阵
%%
if size(y,2)>1 %变成列向量
    y = y';
end

dim = auxdata.dim;
mu = auxdata.mu;

r1=sqrt((y(1)+mu)^2+y(2)^2+y(3)^2);
r2=sqrt((y(1)-1+mu)^2+y(2)^2+y(3)^2);


dydt=[y(4);
    y(5);
    y(6);
    y(1)*(1-(1-mu)/r1^3-mu/r2^3)-mu*(1-mu)*(1/r1^3-1/r2^3)+2*y(5);
    y(2)*(1-(1-mu)/(r1^3)-mu/(r2^3))-2*y(4);
    -y(3)*((1-mu)/(r1^3)+mu/(r2^3))];

if length(y)>dim
    phi = reshape(y(dim+1:end), dim, dim);
    Uxx=1-(1-mu)*(1/r1^3-3*((y(1)+mu)^2)/r1^5)-mu*(1/r2^3-3*((y(1)-1+mu)^2)/r2^5);
    Uyy=1-(1-mu)*(1/r1^3-3*y(2)^2/r1^5)       -mu*(1/r2^3-3*y(2)^2/r2^5);
    Uzz= -(1-mu)*(1/r1^3-3*y(3)^2/r1^5)       -mu*(1/r2^3-3*y(3)^2/r2^5);
    Uxy=3*(1-mu)*(y(1)+mu)*y(2)/r1^5        +3*mu*(y(1)-1+mu)*y(2)/r2^5;
    Uxz=3*(1-mu)*(y(1)+mu)*y(3)/r1^5        +3*mu*(y(1)-1+mu)*y(3)/r2^5;
    Uyz=3*(1-mu)*y(2)*y(3)/r1^5             +3*mu*y(2)*y(3)/r2^5;
    
    UXX=[Uxx,Uxy,Uxz
        Uxy,Uyy,Uyz
        Uxz,Uyz,Uzz];
    
    OMG=[0 2 0
        -2 0 0
        0 0 0];
    
    F=[zeros(3,3),eye(3);
        UXX, OMG];
    
    Dphi=F*phi;
    dydt = [dydt;
        Dphi(:)];
end
