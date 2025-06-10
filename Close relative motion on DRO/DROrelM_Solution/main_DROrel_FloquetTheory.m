%% 基于Floquet theory计算DRO相对运动半解析解
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../subF_eom(CR3BP)')

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


isplot = 1;

%% 积分与计算
% 绝对运动初值
% 初值由共同质心旋转坐标系S下的动力学（陈冠华程序）导出：原点位于共同质心，x轴从地球指向月球，y轴指向月球绕地球旋转方向
% 绝对动力学位于月球旋转坐标系M：原点位于月球，x轴从地球指向月球，y轴指向月球绕地球旋转的方向
% 惯性坐标系：原点位于共同质心，x、y、z轴与初始时刻的M坐标系相同
% 从S至M的位置转换：rM = rS + [1-mu;0]; 速度转换： vM = vS;
% load('InitialState_DRO_2to1')
% load('DRO_all.mat')
% load('DRO_all_SEsyPeriod')
load('DRO_all_MoonPeriod')

numDRO = 4;
x0_DRO_S_2d = state_ini_all(numDRO,:);
x0_DRO_M_3d = [(x0_DRO_S_2d(1:2)-[1-con.mu,0]),0,x0_DRO_S_2d(3:4),0];
para.T0 = J_period_all(numDRO,2)*2*pi; % 轨道周期


%% 线性化相对运动的周期状态转移矩阵（φ(0)^-1 φ(T)）
sol_temp = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[0 1*para.T0], [x0_DRO_M_3d, zeros(1,6), reshape(eye(6),1,36)], opts_ode);
M_REL_lin = reshape(sol_temp.y(13:end,end),6,6);%状态转转移矩阵
[Meigve, Meigva] = eig(M_REL_lin);%特征值
% remove small numbers in Meigve and obtain the accurate eigen vector
Meigve_imga = imag(Meigve); Meigve_imga(abs(Meigve_imga)<1e-5) = 0;
Meigve_real = real(Meigve); Meigve_real(abs(Meigve_real)<1e-5) = 0;
Meigve = Meigve_real+Meigve_imga*1i;

% 判断单位特征值在第1-2行，还是第3-4行
[a,sort_index] = sort(real(diag(Meigva(1:4,1:4))));
index_temp1 = 3*(sort_index(1)>2) + 1*(sort_index(1)<=2);
index_temp2 = 1*(sort_index(1)>2) + 3*(sort_index(1)<=2);
H = blkdiag([real(Meigva(index_temp1,index_temp1)),-imag(Meigva(index_temp1,index_temp1)); ...
    imag(Meigva(index_temp1,index_temp1)),real(Meigva(index_temp1,index_temp1))],...
    [1,1; 0,1],...
    [real(Meigva(5,5)),-imag(Meigva(5,5)); imag(Meigva(5,5)),real(Meigva(5,5))]);
p2_real = real(Meigve(:,index_temp1)); p2_imag = imag(Meigve(:,index_temp1+1));
p3 = real(Meigve(:,index_temp2)); p3 = sign(p3(2)).*p3;
p4 = pinv(M_REL_lin-eye(6))*p3; % 由于M_REL_lin-eye(6)奇异，采用广义逆解
try
    p4 = p4-p4(2)/p3(2)*p3;
catch
end
p3 = p3/norm(p3); 
if norm(p4)>10 % 精度太高有时候算出来的p4超级大
    p4 = pinv(M_REL_lin-eye(6),1e-10)*p3;
end
p4 = p4/norm(p3);

p6_real = real(Meigve(:,6)); p6_imag = imag(Meigve(:,6));

S = [p2_real,p2_imag,p3,p4,p6_real,p6_imag];%特征向量
flag_Jordan = norm(S*H*S^-1-M_REL_lin);

J = 1/para.T0*logm(H);
J = J.*(abs(J)>1e-8);
omega1 = J(2,1);
omega2 = J(3,4); % 正好是1/T0
omega3 = J(6,5);

expJt = @(t)blkdiag([cos(omega1*t),-sin(omega1*t);sin(omega1*t),cos(omega1*t)],...
    [1,omega2*t;0,1],...
    [cos(omega3*t),-sin(omega3*t);sin(omega3*t),cos(omega3*t)]);
flag_Floquet = norm(M_REL_lin * S * expJt(-para.T0) - S);

% resort data
Meigva_diag = diag(Meigva)';
Meigva_diag = Meigva_diag([index_temp1,index_temp1+1,index_temp2,index_temp2+1,5,6]);%特征向量
Meigve = Meigve(:,[index_temp1,index_temp1+1,index_temp2,index_temp2+1,5,6]);

% alpha其实就是算拟周期的时候的旋转量，计算的时候取正虚部的特征值计算
alpha = [atan2(imag(Meigva_diag(2)),real(Meigva_diag(2))),...
        atan2(imag(Meigva_diag(6)),real(Meigva_diag(6)))];
T12_propotion = 2*pi./alpha;

% save data
Sol_linear.vec1 = p2_real;
Sol_linear.vec2 = p2_imag;
Sol_linear.vec3 = p3;
Sol_linear.vec4 = p4;
Sol_linear.vec5 = p6_real;
Sol_linear.vec6 = p6_imag;

% 这里采用第二个和第六个特征值，因为其虚部是正的，如此可以确保alpha大于0
EigenVector.p1_real = p2_real;
EigenVector.p1_imag = p2_imag;
EigenVector.p3 = p3;
EigenVector.p4 = p4;
EigenVector.p5_real = p6_real;
EigenVector.p5_imag = p6_imag;
save ../subF_eom(CR3BP)/FloquetEig12 Meigve Meigva_diag EigenVector Sol_linear state_ini_all J_period_all x0_DRO_M_3d para con
save('Sol_linear.mat','Sol_linear');
save('Meigva_diag.mat','Meigva_diag');
%% 画图
if isplot == 1
    dt = 5*para.T0;
    % dt = 1*T0;
    length_t = 2000;
    t_sample = linspace(0,dt,length_t);
    opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

    for jj_index = 1:6
        eval(['x0_REL = Sol_linear.vec',num2str(jj_index),'*1e-5;'])
        x0_REL =x0_REL';
        % 标称轨道及相对运动 的积分
        sol = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [x0_DRO_M_3d, x0_REL], opts);
        sol_sample = deval(sol,t_sample);
        abs_motion_L_M = sol_sample(1:6,:);
        rel_motion_L_linear = [sol_sample(7:9,:)*con.r_norma;sol_sample(10:12,:)*con.v_norma];

        % 将绝对运动轨道转移至惯性系中
        abs_motion_M = sol_sample(1:6,:);

%         rel_motion_M_linear = T_TCO2TCR_CR3BP(rel_motion_L_linear',abs_motion_M','LVLH',con.mu)';
%         rel_motion_VVLH_linear = T_TCR2TCO_CR3BP(rel_motion_M_linear',abs_motion_M','VVLH',con.mu)';

        figure(1)
        subplot(3,2,jj_index)
        if jj_index<5
            plot(rel_motion_L_linear(1,1),rel_motion_L_linear(2,1),'g^'); hold on
            plot(rel_motion_L_linear(1,end),rel_motion_L_linear(2,end),'r*'); 
            plot(rel_motion_L_linear(1,:),rel_motion_L_linear(2,:),'Color',[0, 114, 189]/255,'LineWidth',1.5)
            hold off
    %         legend('InitialPos','FinalPos')
            legend('初始时刻','终端时刻')
            xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]')
            axis equal;
            x_ratio = 3.3; % x轴显示与图形中点的比例
            y2x_ratio = 0.7; % 画图时y轴显示与x轴显示的比例
            x_max = max(rel_motion_L_linear(1,:));
            x_min = min(rel_motion_L_linear(1,:));
            x_middle = (x_max+x_min)/2;
            x_diff = x_max - x_min;
            y_middle = (max(rel_motion_L_linear(2,:))+min(rel_motion_L_linear(2,:)))/2;
            xlim([x_middle - x_ratio/2*x_diff, x_middle + x_ratio/2*x_diff]) 
            ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 
        else
            plot(t_sample/para.T0,rel_motion_L_linear(3,:),'LineWidth',1.5)
            xlim([min(t_sample/para.T0), max(t_sample/para.T0)]) 
            xlabel('\itt \rm[T]'); ylabel('\itz \rm[km]');
        end
           
        grid on; grid minor

    end
end