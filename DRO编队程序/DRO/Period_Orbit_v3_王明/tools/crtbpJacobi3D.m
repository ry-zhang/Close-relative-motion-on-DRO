
function CC = crtbpJacobi3D (xx, mu)
% 计算雅可比能量C
%
% 输入：
% xx       [nx4]
% 输出：
% C        [nx1]
%
% 2015/4/10, v0
% 2020/02/05, v1：xx: n * 6，匹配ode113
% Copyright(C) Chen Zhang
% -------------------------------------------------------------
if size(xx , 2) ~= 6
    disp('xx是n * 6'); 
end
r1 = sqrt((xx(: , 1) + mu).^2   + xx(: , 2).^2 + xx(: , 3).^2) ;
r2 = sqrt((xx(: , 1) + mu - 1).^2 + xx(: , 2).^2 + xx(: , 3).^2) ;
CC = -(xx(: , 4).^2 + xx(: , 5).^2 + xx(: , 6).^2) ...
    + (xx(: , 1).^2 + xx(: , 2).^2) ...
    + 2 * (1 - mu) ./ r1 + 2 * mu ./ r2;
end
