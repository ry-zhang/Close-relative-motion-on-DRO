%% 分析不同阶次的傅里叶级数的计算量和精度
% 2022-12-3
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
%% 
clear
addpath('../../subF_eom(CR3BP)')

format longg
format compact

set(0,'defaultAxesFontName', 'TimesSimSun','defaultTextFontName', 'TimesSimSun');
set(0,'defaultAxesFontSize',15,'defaultTextFontSize',15)
set(0,'defaultLineLineWidth',1.5)

%% 常数与变量
load('FloquetEig12.mat')
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);
para.x0_DRO = x0_DRO_M_3d;

%% 相对运动积分

dt = 1*para.T0; % 积分时间
length_t = 2000;
t_sample = linspace(0,dt,length_t);
t_sample_day = t_sample*con.T_norma_day;

% 周期轨道
% 标称轨道及相对运动 p3 的积分
sol3 = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [para.x0_DRO, EigenVector.p3'], opts);
sol3_sample = deval(sol3,t_sample);
e3_hat = sol3_sample(7:12,:);

% 平面拟周期
% 这里采用第二个特征值，因为其虚部是正的，如此可以确保alpha大于0
alpha1 = atan2(imag(Meigva_diag(2)),real(Meigva_diag(2)));
para.T1 = 2*pi*para.T0/alpha1;
% 标称轨道及相对运动 p1_real 的积分
sol1 = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [para.x0_DRO, EigenVector.p1_real'], opts);
sol1_sample = deval(sol1,t_sample);
abs_motion_M = sol1_sample(1:6,:);
rel_motion_L_linear_1 = sol1_sample(7:12,:);
% 标称轨道及相对运动 p1_imag 的积分
sol2 = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [para.x0_DRO, EigenVector.p1_imag'], opts);
sol2_sample = deval(sol2,t_sample);
rel_motion_L_linear_2 = sol2_sample(7:12,:);

% e_1
e1_hat = cos(-alpha1*t_sample/para.T0).*rel_motion_L_linear_1 + sin(-alpha1*t_sample/para.T0).*rel_motion_L_linear_2;
e2_hat = -sin(-alpha1*t_sample/para.T0).*rel_motion_L_linear_1 + cos(-alpha1*t_sample/para.T0).*rel_motion_L_linear_2;

% 法向拟周期
% 这里采用第六个特征值，因为其虚部是正的，如此可以确保alpha大于0
alpha2 = atan2(imag(Meigva_diag(6)),real(Meigva_diag(6)));
para.T2 = 2*pi*para.T0/alpha2;
% 标称轨道及相对运动 p5_real 的积分
sol5 = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [para.x0_DRO, EigenVector.p5_real'], opts);
sol5_sample = deval(sol5,t_sample);
rel_motion_L_linear_5 = sol5_sample(7:12,:);
% 标称轨道及相对运动 p5_imag 的积分
sol6 = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 dt], [para.x0_DRO, EigenVector.p5_imag'], opts);
sol6_sample = deval(sol6,t_sample);
rel_motion_L_linear_6 = sol6_sample(7:12,:);

% e_5
e5_hat = cos(-alpha2*t_sample/para.T0).*rel_motion_L_linear_5 + sin(-alpha2*t_sample/para.T0).*rel_motion_L_linear_6;
e6_hat = -sin(-alpha2*t_sample/para.T0).*rel_motion_L_linear_5 + cos(-alpha2*t_sample/para.T0).*rel_motion_L_linear_6;

%% 用傅里叶变换拟合平面内拟周期
% e1_hat、e2_hat
N_all = 2:2:100;
t_ode_all = zeros(size(N_all));
t_FFT_all = t_ode_all;
error_max_rel_all = zeros(2,length(N_all));

for jj_loop = 1:length(N_all)
    N = N_all(jj_loop);
    disp(N); 
    coe.N = N;
    t_sample_fit = linspace(0,para.T0*(coe.N-1)/coe.N,coe.N);

    % %% 用傅里叶变换拟合平面内周期轨道
    % t_sample2 = rand(10,1);
    % e3_hat
    tic
    for ii = 1:10
        sol3_sample = deval(sol3,t_sample);
    end
    t_ode_all(jj_loop) = toc/10;

    sol3_sample_fit = deval(sol3,t_sample_fit);
    e3_hat_fitsam = sol3_sample_fit(7:12,:)';


    coe.c1_e3hat = DFTmatrix(coe.N)*e3_hat_fitsam;

    theta0_all_fit = linspace(0,2*pi,length_t);
    % theta0_all_fit = t_sample_fit*2*pi;
    tic
    e3_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e3hat);
    for ii = 1:100
        e3_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e3hat);
    end
    t_FFT_all(jj_loop) = toc/100;

    error3 = e3_refit-e3_hat';
    error3_max = max(abs(error3(:,[1,2,4,5])));
    error3_max_rel = max(abs(error3(:,[1,2,4,5])))./mean(abs(e3_hat([1,2,4,5],:)'));
    error_max_rel_all(:,jj_loop) = error3_max_rel(1:2);
    % disp('error3:'); 
    % 绝对误差
    % disp(max(abs(error3(:,[1,2,4,5]))))
    % 相对误差
    % disp(max(error3_rel.*(abs(e3_hat([1,2,4,5],:)')>1e-3)))
    % disp(max(abs(error3(:,[1,2,4,5])))./mean(abs(e3_hat([1,2,4,5],:)')))
end

% 画图
order_all = N_all/2;
figure(1)
plot(order_all,t_FFT_all)
xlabel('截断阶次'); ylabel('计算时间[sec]')
exportgraphics(gcf,'计算时间.jpg','Resolution',600)
figure(2)
semilogy(order_all,error_max_rel_all)
xlabel('截断阶次'); ylabel('最大相对精度')
legend('\itx','\ity')
exportgraphics(gcf,'计算精度.jpg','Resolution',600)