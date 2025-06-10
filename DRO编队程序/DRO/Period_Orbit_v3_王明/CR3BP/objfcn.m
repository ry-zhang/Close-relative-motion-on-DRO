function [fobj , Gobj] = objfcn(yy , aux)
%
% Objective function that considers the final point of the trajectory
% (xxn) as a NLP variable (i.e., there is no final integration to impose
% the right boundary condition).
%
% 作者：张晨
% 2021年5月30日
%%%%%%%%%%%%%%%%%%%%%%%%

% number of patch points
n = aux.node_n;

% objective function
fobj = 1;

if nargout > 1 % compute gradient if necessary
    
    Gobj = zeros(1 , 6 * n + 1);
    
end

end