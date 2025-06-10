
function plot_solution(xj , Gc_eq , aux)
% 
% 画仿真结果
%
% 作者：张晨
% 2021年6月8日
%%%%%%%%%%%%%%%%%%%%%%%%

% mu = aux.ems_mu;
 LU = aux.LU;
% VU = aux.ems_VU;
% TU = aux.ems_TU;


if aux.UseSymmetric_IO == 1
    resoOrb_xi = [xj(1), 0, xj(2), 0, xj(3), 0]';
    resoOrb_P = 2*xj(end);
else
    resoOrb_xi = xj(1:6);
    resoOrb_P = xj(end);
end
% ------------------------- rot frame -------------------------
h2 = figure(2); hold on; grid on;
set(h2 , 'position' , [750 , 100 , 600 , 500]);
crtbpMarkEM;
[tt_rot , xx_rot] = crtbpFlow3D(resoOrb_xi , resoOrb_P , 'b' , 1 , aux);
tt_interp = linspace(0 , resoOrb_P , aux.node_n + 1)';
xx_rot_interp = interp1(tt_rot , xx_rot , tt_interp' , 'spline') ;
plot3(xx_rot_interp(1:end - 1 , 1) , xx_rot_interp(1:end - 1 , 2) , xx_rot_interp(1:end - 1 , 3) , 'bx' , 'linewidth' , 1);
plot3(xx_rot_interp(1 , 1) , xx_rot_interp(1 , 2) , xx_rot_interp(1 , 3) , 'rx' , 'linewidth' , 2 , 'markersize' , 8);
axis equal;



% ------------------------- ini frame -------------------------
h1 = figure(1); hold on; grid on;
set(h1 , 'position' , [100 , 100 , 600 , 500]);

plot_o([0 , 0 , 0] , 1 , 'k--');
text(0,0,0 , 'Earth');
text(1,0,0 , 'Moon');

[xs , ys , zs] = sphere(30) ;
surf(6378 * xs / LU,  6378 * ys / LU, 6378 * zs / LU, 'FaceColor', 'w') ; hold on;
surf(1738 * xs / LU + 1 , 1738 * ys / LU, 1738 * zs / LU, 'FaceColor', 'w') ; hold on;

[xx_ini , tt_ini] = crtbpSynodic2P1centered3D(xx_rot , tt_rot , aux.mu);
plot3(xx_ini(: , 1) , xx_ini(: , 2) , xx_ini(: , 3) , 'b' , 'linewidth' , 1);
tt_interp = linspace(0 , resoOrb_P , aux.node_n + 1)';
xx_ini_interp = interp1(tt_ini , xx_ini , tt_interp' , 'spline') ;
plot3(xx_ini_interp(1:end - 1 , 1) , xx_ini_interp(1:end - 1 , 2) , xx_ini_interp(1:end - 1 , 3) , 'bx' , 'linewidth' , 1);
plot3(xx_ini_interp(1 , 1) , xx_ini_interp(1 , 2) , xx_ini_interp(1 , 3) , 'bx' , 'linewidth' , 2 , 'markersize' , 8);

axis equal;


% ------------------------- jacobian structure -------------------------
% h3 = figure(3); hold on; grid on;
% set(h3 , 'position' , [300 , 300 , 600 , 500]);
% spy(Gc_eq)


% ------------------------- disp -------------------------
fprintf(' --------------------------------------------------- \n')
fprintf('周期轨道初值：\n')
resoOrb_xi

fprintf('周期轨道周期：\n')
resoOrb_P


end
