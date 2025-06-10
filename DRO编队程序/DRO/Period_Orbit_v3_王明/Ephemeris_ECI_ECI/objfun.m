function  [fobj , Gobj] = objfun(yy , aux)
% number of vars
nVar = size(yy , 1);
N = aux.node_n;
% objective function
fobj = norm(yy(1:6) - yy(6*(N-1)+1:6*(N-1)+6));

if nargout > 1 % compute gradient if necessary
    
    Gobj = zeros(1 , nVar);
    Gobj (1 , 1: 6) = yy(1:6)'/fobj ;
    Gobj (1 , 6*(N-1)+1:6*(N-1)+6) = -yy(6*(N-1)+1:6*(N-1)+6)'/fobj ;
end
end

