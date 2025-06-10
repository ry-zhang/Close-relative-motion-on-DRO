set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'DefaultAxesFontsize',15);
set(0,'DefaultTextFontsize',15);
set(0,'defaultLineLineWidth',1.5)

load('FloquetEig12')

r_normal_all = [0,0,1; 0,1,0; 1,0,0; 0,1,1];
r_cent_all = [0,0,0; 0,0,1; 0,1,0; 1,0,0];
% r_cent_all = [1,0,0; 2,0,0; 5,0,0; 10,0,0];
% r_cent_all = [0,0,0; 3*[0.707,0.707,0]; 7*[0.707,0.707,0]; 12*[0.707,0.707,0]];
radi_all = [1,3,10];

% flag_str = '法向量'; size_r = size(r_normal_all,1);
flag_str = '圆心'; size_r = size(r_cent_all,1);
% flag_str = '半径'; size_r = size(radi_all,1);
for ii_loop = 1:size_r
    n_samp = 800; % 每圈的点数
    dt = 10/24/con.T_norma_day/(n_samp-1); % 转移时间
    % 空间圆方程
    if strcmp(flag_str,'法向量')
        radi = 1; % 半径，km
        r_cent = [0,0,0]; % 圆心，km
        r_normal = r_normal_all(ii_loop,:); % 法向量
    elseif strcmp(flag_str,'圆心')
        radi = 1; % 半径，km
        r_cent = r_cent_all(ii_loop,:); % 法向量
        r_normal = [0,1,1]; % 法向量
    elseif strcmp(flag_str,'半径')
        radi = radi_all(ii_loop); % 半径，km
        r_cent = [0,0,0]; % 法向量
        r_normal = [0,1,1]; % 法向量
    end
    r_normal = r_normal/norm(r_normal); % 单位法向量
    % 计算r_a与r_b,先对r_normal升序排列，以寻找最大的非零值，预防奇异
    [r_normal_resort,r_normal_index] = sort(abs(r_normal),'ascend');
    r_a1 = 0.1;
    r_a2 = 0.2;
    r_a3 = -(r_a1*r_normal_resort(1)+r_a2*r_normal_resort(2))/r_normal_resort(3);
    r_a_temp = [r_a1,r_a2,r_a3];
    eye3 = eye(3);
    r_a = r_a_temp*eye3(r_normal_index,:)/norm(r_a_temp);
    r_b = cross(r_normal,r_a);
    theta_all = linspace(0,2*pi*n_samp/(n_samp+1),n_samp);
    r_chaser_km = r_cent + radi*cos(theta_all)'*r_a + radi*sin(theta_all)'*r_b;
    
    plot3(r_chaser_km(:,1),r_chaser_km(:,2),r_chaser_km(:,3),'LineWidth',2); hold on;
%     plot3(r_chaser_km(:,1),r_chaser_km(:,2),r_chaser_km(:,3),'Color',[200,200,200]/255,'LineWidth',2); hold on;
end
hold off

xlabel('{\itx_L} [km]'); ylabel('{\ity_L} [km]'); zlabel('{\itz_L} [km]');
axis equal; 
% title('参考空间圆轨道'); 
grid on; grid minor
box on
% set(gcf,'Renderer', 'painters')

if strcmp(flag_str,'法向量')
    r_all = r_normal_all;
elseif strcmp(flag_str,'圆心')
    r_all = r_cent_all;
elseif strcmp(flag_str,'半径')
    r_all = radi_all;
end
str = strings(1,size(r_all,1));
for jj_loop = 1:size(r_cent_all,1)
    str(jj_loop) = [flag_str,' = [',num2str(r_all(jj_loop,1)),',',...
        num2str(r_all(jj_loop,2)),',',num2str(r_all(jj_loop,3)),']^T'];
end
legend(str,'Location','eastoutside')

% set(gcf,'Color',[255,255,255]/255);
% export_fig RefOrbit.png -r600