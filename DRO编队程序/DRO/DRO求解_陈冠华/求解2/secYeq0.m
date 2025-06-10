function [value,isterm,drct] = secYeq0(t,X)
% event defined by y's value

value = 1;
secY = 0; 
for ii = 1:length(secY)
    value = value*(X(2)-secY(ii));
end

isterm = 0;
drct = 1;
