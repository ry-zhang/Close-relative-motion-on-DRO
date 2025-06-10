clear
addpath('../../subF_eom(CR3BP)')
load('generalSolFFT_12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 15);%坐标轴
set(0,'defaultTextFontSize', 15);%文字
set(0,'defaultLineLineWidth',1.5)

%% 计算ECI下的相对运动
dt = 7*para.T0; % 积分时间                  
length_t = 500;
t_sample = linspace(0,dt,length_t);
% t_sample = linspace(para.T0,dt,length_t);
ratio = 0.9/con.r_norma;

% 自然轨迹 绝对运动
sol_DRO = ode113(@(t,x)eom_abs3b(t,x,con.mu),[0 dt], para.x0_DRO, opts);
x_chief_MCR = deval(sol_DRO,t_sample);
% 将月球为中心天体的旋转坐标系转化到地球为中心天体的旋转坐标系，绝对运动转换
x_chief_ECR = x_chief_MCR; x_chief_ECR(1:3,:) = x_chief_MCR(1:3,:)-[-1,0,0]';
x_chief_ECI = synodic2inertial(x_chief_ECR,t_sample);
theta1_0_all = linspace(0,2*pi,100); % 0.625*pi

for ii_loop = 1:3
    switch ii_loop
        case 1 
            theta1_0 = 1.01546429206943;
            k1 = 0; k2 = 0; k0 = 1; theta2_0 = 1.266*pi;
        case 2
            theta1_0 = 1.01546429206943;
            k1 = 1.3; k2 = 0; k0 = 0; theta2_0 = 1.266*pi;
        case 3
            theta1_0 = 1.01546429206943;
            k1 = 0; k2 = 1; k0 = 0; theta2_0 = 1.266*pi;
    end
    
    x_rel_L = generalSol_relMotion(t_sample,k0,k1,k2,theta1_0,theta2_0,para,coe);
    x_rel_MCR = T_TCO2TCR_CR3BP(x_rel_L',x_chief_MCR','LVLH',con.mu)';

    % 放缩
    r_rel_L = ratio*con.r_norma*x_rel_L(1:3,:);
    
    % 地球为中心天体的惯性坐标系 下的轨道坐标系
    % natural 轨道段
    x_rel_ECR = x_rel_MCR; % 旋转坐标系下的相对运动是相同的，与旋转坐标系质心无关
    x_rel_ECI = synodic2inertial(x_rel_ECR,t_sample);
    x_rel_ECIL = T_TCI2TCO_E_CR3BP(x_rel_ECI',x_chief_ECI',t_sample,'LVLH',con.mu)';
    % 由于DRO绕月逆行，因此LVLH的z轴方向与MCR的z轴相反。
    
    r_rel_ECIL_km = ratio*con.r_norma*x_rel_ECIL(1:3,:);
    
    if ii_loop < 3
        % 画MCR LVLH中的编队
        figure(1)
        plot3(r_rel_L(1,:),r_rel_L(2,:),r_rel_L(3,:),'LineWidth',1.5);
        hold on
        plot3(r_rel_L(1,1),r_rel_L(2,1),r_rel_L(3,1),'g^');
        plot3(r_rel_L(1,end),r_rel_L(2,end),r_rel_L(3,end),'rv');
        plot(0,0,'ks','MarkerSize',5)
        xlabel('\itx_L \rm'); ylabel('\ity_L \rm'); 
        hold off; axis equal; grid on; box on
        xlim([-0.401,0.401]); ylim([-0.601,0.601]); zlim([-0.801,0.801]);
        if ii_loop == 1
            title('平面周期模态')
        elseif ii_loop == 2
            title('平面拟周期模态')
        end
        view([0,90])
        export_fig LVLHMCR.png -r600
        
        % 画ECI LVLH中的编队
        figure(2)
        plot3(r_rel_ECIL_km(1,:),r_rel_ECIL_km(2,:),r_rel_ECIL_km(3,:),'LineWidth',1.5);
        hold on
        p1 = plot3(r_rel_ECIL_km(1,1),r_rel_ECIL_km(2,1),r_rel_ECIL_km(3,1),'g^');
        p2 = plot3(r_rel_ECIL_km(1,end),r_rel_ECIL_km(2,end),r_rel_ECIL_km(3,end),'rv');
        plot(0,0,'ks','MarkerSize',5) 
        xlabel('\itx_{LE}'); ylabel('\ity_{LE}');
        hold off; axis equal; grid on; box on
        xlim([-0.401,0.401]); ylim([-0.601,0.601]); zlim([-0.801,0.801]);
        if ii_loop == 1
            title('平面周期模态')
        elseif ii_loop == 2
            title('平面拟周期模态')
        end
        view([0,90])
        export_fig LVLHECI.png -r600
    else
        % 画MCR LVLH中的编队
        figure(1)
        plot(t_sample,r_rel_L(3,:),'LineWidth',1.5);
        xlabel('{\itt} {[{\itT}_0]}'); ylabel('\itz_{L}');
        title('法向拟周期模态'); grid on;
        export_fig LVLHMCR.png -r600
        % 画ECI LVLH中的编队
        figure(2)
        plot(t_sample,r_rel_ECIL_km(3,:),'LineWidth',1.5);
        xlabel('{\itt} {[{\itT}_0]}'); ylabel('\itz_{LE}');
        title('法向拟周期模态'); grid on;
        export_fig LVLHECI.png -r600
    end
end
