function [value,isterm,drct] = secYeq0Stop(t,X)

value = 1;
secY = 0; % ps defined here as y = 0
for ii = 1:length(secY)
    value = value*(X(2)-secY(ii));
end

isterm = 1;
drct = 1;
