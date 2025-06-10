%% 计算周期状态转移矩阵
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../../subF_eom(CR3BP)')

%% 常数与变量
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2
% con.mu = 0.01211; % 
con.mu = 0.01215; % 20200531
con.r_norma = 3.84399*10^5; % km
% T_M = 27.321661; % day
con.T_norma = sqrt(con.r_norma^3/(mu_E+mu_M)); % s
con.T_norma_day = con.T_norma/3600/24;
con.v_norma = con.r_norma/con.T_norma; % km/s
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);
plot_flag = 0;

%% 积分与计算
% 绝对运动初值
% 初值由共同质心旋转坐标系S下的动力学（陈冠华程序）导出：原点位于共同质心，x轴从地球指向月球，y轴指向月球绕地球旋转方向
% 绝对动力学位于月球旋转坐标系M：原点位于月球，x轴从地球指向月球，y轴指向月球绕地球旋转的方向
% 惯性坐标系：原点位于共同质心，x、y、z轴与初始时刻的M坐标系相同
% 从S至M的位置转换：rM = rS + [1-mu;0]; 速度转换： vM = vS;
% load('DRO_3to1')
% load('DRO_all3.mat')
% load('DRO_all4.mat')
% flag = '_resonance';
flag = '4-2';
load(['DRO_all',flag,'.mat'])

num = size(state_ini_all,1);

flag_Jordan_all = zeros(num,1);
flag_Floquet_all = zeros(num,1);
Meigve_all = zeros(6,6,num);
Meigva_all = zeros(num,6);
alpha_all = zeros(num,2);
EigenVector_all(num).p1_real = [];
x0_DRO_M_3d_all = zeros(num,6);
DRO_error = zeros(num,1);
RELperi_error = zeros(num,1);
for i_index = 1:num
% for i_index = 96:num
    disp(num2str(i_index))
    x0_DRO_S_2d = state_ini_all(i_index,:);
    x0_DRO_M_3d = [(x0_DRO_S_2d(1:2)-[1-con.mu,0]),0,x0_DRO_S_2d(3:4),0];
    x0_DRO_M_3d_all(i_index,:) = x0_DRO_M_3d;
    T0 = J_period_all(i_index,2)*2*pi; % 轨道周期
    %% 线性化相对运动的状态转移矩阵（φ(0)^-1 φ(T)）
    sol_temp = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[0 1*T0], [x0_DRO_M_3d, zeros(1,6), reshape(eye(6),1,36)], opts_ode);
    M_REL_lin = reshape(sol_temp.y(13:end,end),6,6);
    [Meigve, Meigva] = eig(M_REL_lin);
    % remove small numbers in Meigve and obtain the accurate eigen vector
    Meigve_imga = imag(Meigve); Meigve_imga(abs(Meigve_imga)<1e-5) = 0;
    Meigve_real = real(Meigve); Meigve_real(abs(Meigve_real)<1e-5) = 0;
    Meigve = Meigve_real+Meigve_imga*1i;
    if imag(Meigva(5,5)) == 0
        [~,ve_sort] = sort(real(diag(Meigva(5:6,5:6))));
        Meigva(5:6,5:6) = Meigva(4+ve_sort,4+ve_sort);
        Meigve(:,5:6) = Meigve(:,4+ve_sort);
    end
    
    % 判断单位特征值在第1-2行，还是第3-4行
    [a,sort_index] = sort(real(diag(Meigva(1:4,1:4))));
    index_temp1 = 3*(sort_index(1)>2) + 1*(sort_index(1)<=2);
    index_temp2 = 1*(sort_index(1)>2) + 3*(sort_index(1)<=2);
    H = blkdiag([real(Meigva(index_temp1,index_temp1)),-imag(Meigva(index_temp1,index_temp1)); ...
        imag(Meigva(index_temp1,index_temp1)),real(Meigva(index_temp1,index_temp1))],...
        [1,1; 0,1],...
        [real(Meigva(5,5)),-imag(Meigva(5,5)); imag(Meigva(5,5)),real(Meigva(5,5))]);
    p2_real = real(Meigve(:,index_temp1)); p2_imag = imag(Meigve(:,index_temp1+1));
    if p2_real(2) ~= 0
        p2_imag = p2_imag*sign(p2_real(2)); p2_real = p2_real*sign(p2_real(2));
    else
        % 注意，这种情况下，实际上将theta2平移了pi/2
        p2_imag = real(Meigve(:,index_temp1)); p2_real = imag(Meigve(:,index_temp1));
    end
    p3 = real(Meigve(:,index_temp2)); p3 = sign(p3(2)).*p3; % y轴正向
    p4 = pinv(M_REL_lin-eye(6))*p3; % 由于M_REL_lin-eye(6)奇异，采用广义逆解
    p3 = p3/norm(p3); 
    if norm(p4)>10 % 精度太高有时候算出来的p4超级大
        p4 = pinv(M_REL_lin-eye(6),1e-10)*p3;
    end
    p4 = p4/norm(p3);

    p6_real = real(Meigve(:,6)); p6_imag = imag(Meigve(:,6));
    if p6_real(6) ~= 0
        p6_imag = p6_imag*sign(p6_real(6)); p6_real = p6_real*sign(p6_real(6)); 
    end
    S = [p2_real,p2_imag,p3,p4,p6_real,p6_imag];
    flag_Jordan = norm(S*H*S^-1-M_REL_lin);

    J = 1/T0*logm(H);
    J = J.*(abs(J)>1e-8);
    omega1 = J(2,1);
    omega2 = J(3,4); % 正好是1/T0
    omega3 = J(6,5);

    expJt = @(t)blkdiag([cos(omega1*t),-sin(omega1*t);sin(omega1*t),cos(omega1*t)],...
        [1,omega2*t;0,1],...
        [cos(omega3*t),-sin(omega3*t);sin(omega3*t),cos(omega3*t)]);
    flag_Floquet = norm(M_REL_lin * S * expJt(-T0) - S);
    
    % save data
    flag_Jordan_all(i_index) = flag_Jordan;
    flag_Floquet_all(i_index) = flag_Floquet;
    Meigva_diag = diag(Meigva)';
    Meigva_all(i_index,:) = Meigva_diag([index_temp1,index_temp1+1,index_temp2,index_temp2+1,5,6]);
    Meigve_all(:,:,i_index) = S;
    EigenVector_all(i_index).p1_real = p2_real;
    EigenVector_all(i_index).p1_imag = p2_imag;
    EigenVector_all(i_index).p3 = p3;
    EigenVector_all(i_index).p4 = p4;
    EigenVector_all(i_index).p5_real = p6_real;
    EigenVector_all(i_index).p5_imag = p6_imag;
    % alpha其实就是算拟周期的时候的旋转量，计算的时候取正虚部的特征值计算
    alpha_all(i_index,:) = [atan2(imag(Meigva_all(i_index,2)),real(Meigva_all(i_index,2))),...
        atan2(imag(Meigva_all(i_index,6)),real(Meigva_all(i_index,6)))];

    %% 画图
    if plot_flag == 1
        dt = 10*T0; % 积分时间
        % dt = 1*T0; % 积分时间
        length_t = 2000;
        t_sample = linspace(0,dt,length_t);

        for jj_index = 1:6
            x0_REL = S(:,jj_index) *1e-7;
            x0_REL = x0_REL';
            % 标称轨道及相对运动 的积分
            sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts_ode);
            sol_sample = deval(sol,t_sample);
            abs_motion_L_M = sol_sample(1:6,:);
            rel_motion_L_linear = sol_sample(7:12,:);
            
            figure(1)
            subplot(3,2,jj_index)
            if jj_index<5
                plot(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:),'LineWidth',1.5)
                xlabel('x_L'); ylabel('y_L');
                axis equal;
                x_ratio = 3.3; % x轴显示与图形中点的比例
                y2x_ratio = 0.7; % 画图时y轴显示与x轴显示的比例
                xyLim(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:),x_ratio,y2x_ratio)
            else
                plot(t_sample/T0,rel_motion_L_linear(3,:),'LineWidth',1.5)
                xlabel('t [T]'); ylabel('z_L');
            end
            title(['x_',num2str(jj_index),'(t)']); set(gca,'FontSize',13)
            grid on; grid minor

        end
        
    end
    
    DRO_error(i_index) = norm(sol_temp.y(1:6,1)-sol_temp.y(1:6,end));
    sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 T0], [x0_DRO_M_3d, EigenVector_all(i_index).p3'], opts_ode);
    RELperi_error(i_index) = norm(sol.y(7:12,1)-sol.y(7:12,end));
end

m2n_all = 2*pi./alpha_all;
if plot_flag == 1
    figure(10)
    plot(J_period_all(:,2),m2n_all(:,1),'LineWidth',1.5);
    set(gca,'FontSize',13)
    xlabel('T_{DRO} [TU]'); ylabel('m_1/n_1');
    grid on; grid minor
end

%% 保存数据
% flag_Jordan_all flag_Floquet_all
save(['../../subF_eom(CR3BP)/FloquetEig_all',flag], 'Meigva_all', 'Meigve_all',...
    'EigenVector_all', 'x0_DRO_M_3d_all', 'J_period_all', 'alpha_all', 'm2n_all', 'con')
