
function plot_o(pos , radius , type)
% 
% 画圆圈
% 
% 作者：张晨
% 时间：2021年1月21日
% 单位：中科院空间应用工程与技术中心，空间探索室
% 邮箱：chenzhang.buaa@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%

theta = linspace(0 , 2*pi , 60);
xx = radius * cos(theta) + pos(1);
yy = radius * sin(theta) + pos(2);
zz = 0 * zeros(size(theta)) + pos(3);
plot3(xx , yy , zz , type);

end

