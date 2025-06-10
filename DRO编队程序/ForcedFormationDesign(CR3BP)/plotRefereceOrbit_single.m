set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'DefaultAxesFontsize',15);
set(0,'DefaultTextFontsize',15);
set(0,'defaultLineLineWidth',1.5)

load('FloquetEig12')

        radi = 1; % 半径，km
        r_cent = [0,0,0]; % 法向量
        r_normal = [0.9,0,-0.3]; % 法向量

    n_samp = 800; % 每圈的点数
    dt = 10/24/con.T_norma_day/(n_samp-1); % 转移时间
    % 空间圆方程
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
    p1 = plot3(r_chaser_km(:,1),r_chaser_km(:,2),r_chaser_km(:,3),'k','LineWidth',1.5); hold on;
    r_DRO = [-10,10; 0,0; 0,0]';
    p0 = plot3(r_DRO(:,1),r_DRO(:,2),r_DRO(:,3),'Color',[0, 114, 189]/255,'LineWidth',2); hold off
%     plot3(r_chaser_km(:,1),r_chaser_km(:,2),r_chaser_km(:,3),'Color',[200,200,200]/255,'LineWidth',2); hold on;

axis equal; 
xlim([-2,2]); ylim([-1,1]); zlim([-1,1])
xlabel('{\itx_M} [km]'); ylabel('{\ity_M} [km]'); zlabel('{\itz_M} [km]');
% title('参考空间圆轨道'); 
grid on;
box on
legend([p0,p1],{'DRO','线性化相对运动'},'Location','north')

set(gcf,'Color',[255,255,255]/255);
export_fig RefOrbit.png -r600