%% 几类坐标转换矩阵
% % 从stk导出结果
% icrf_x = [ 1.0000000000     0.0000000001    -0.0000002089 ]  ;
% icrf_z = [ 0.0000002089    0.0000000537    1.0000000000];
% icrf_y = cross(icrf_z, icrf_x);
% rot = [icrf_x;
%     icrf_y
%     icrf_z]';

% 
% ICRF_2_MeanEclpJ2k = [1,-5.09103560100084e-11,2.09100000000000e-07;
%     -8.31000000000000e-08,0.917482040700040,0.397777205100000;
%     -1.91800000000000e-07,-0.397777205100017,0.917482040700000];

obliquity = 23.43662*pi/180;
ICRF_2_MeanEclpJ2k = [1   0   0;
    0    cos(obliquity)   sin(obliquity);
    0   -sin(obliquity)   cos(obliquity)];

% MEME
ICRF_2_J2K = eye(3);


