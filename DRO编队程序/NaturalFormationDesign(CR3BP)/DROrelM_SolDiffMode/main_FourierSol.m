%% CR3BP下DRO的傅里叶形式通解
% 2021-8-30
% by Yang Chihang
% email: ychhtl@foxmail.com
close all
clear
addpath('../../subF_eom(CR3BP)')

format longg
format compact

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
coe.N = 50;
t_sample_fit = linspace(0,para.T0*(coe.N-1)/coe.N,coe.N);

sol1_sample_fit = deval(sol1,t_sample_fit);
sol2_sample_fit = deval(sol2,t_sample_fit); 

e1_hat_fitsam = cos(-alpha1*t_sample_fit/para.T0).*sol1_sample_fit(7:12,:) + sin(-alpha1*t_sample_fit/para.T0).*sol2_sample_fit(7:12,:);
e2_hat_fitsam = -sin(-alpha1*t_sample_fit/para.T0).*sol1_sample_fit(7:12,:) + cos(-alpha1*t_sample_fit/para.T0).*sol2_sample_fit(7:12,:);
e1_hat_fitsam = e1_hat_fitsam'; e2_hat_fitsam = e2_hat_fitsam';

coe.c1_e1hat = DFTmatrix(coe.N)*e1_hat_fitsam;
coe.c1_e2hat = DFTmatrix(coe.N)*e2_hat_fitsam;

theta0_all_fit = linspace(0,2*pi,length_t);
% theta_all_fit = t_sample_fit*2*pi/para.T0;
e1_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e1hat);
e2_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e2hat);

% error
error1 = e1_refit-e1_hat'; 
error2 = e2_refit-e2_hat';
disp('error1:'); disp(max(abs(error1(:,[1,2,4,5]))))
disp('error2:'); disp(max(abs(error2(:,[1,2,4,5]))))

%% 用傅里叶变换拟合平面内周期轨道
% e3_hat
sol3_sample_fit = deval(sol3,t_sample_fit);
e3_hat_fitsam = sol3_sample_fit(7:12,:)';

coe.c1_e3hat = DFTmatrix(coe.N)*e3_hat_fitsam;

theta0_all_fit = linspace(0,2*pi,length_t);
e3_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e3hat);

error3 = e3_refit-e3_hat';
disp('error3:'); 
% 绝对误差
disp(max(abs(error3(:,[1,2,4,5]))))
% 相对误差
disp(max(abs(error3(:,[1,2,4,5])./e3_hat([1,2,4,5],:)').*abs(e3_hat([1,2,4,5],:)'>1e-5)))
%% 用傅里叶变换拟合平面外拟周期
% e5_hat与e6_hat
sol5_sample_fit = deval(sol5,t_sample_fit);
sol6_sample_fit = deval(sol6,t_sample_fit); 

e5_hat_fitsam = cos(-alpha2*t_sample_fit/para.T0).*sol5_sample_fit(7:12,:) + sin(-alpha2*t_sample_fit/para.T0).*sol6_sample_fit(7:12,:);
e6_hat_fitsam = -sin(-alpha2*t_sample_fit/para.T0).*sol5_sample_fit(7:12,:) + cos(-alpha2*t_sample_fit/para.T0).*sol6_sample_fit(7:12,:);
e5_hat_fitsam = e5_hat_fitsam'; e6_hat_fitsam = e6_hat_fitsam';

coe.c1_e5hat = DFTmatrix(coe.N)*e5_hat_fitsam;
coe.c1_e6hat = DFTmatrix(coe.N)*e6_hat_fitsam;

theta0_all_fit = linspace(0,2*pi,length_t);
e5_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e5hat);
e6_refit = real(iDFTmatrix_theta(coe.N,theta0_all_fit) * coe.c1_e6hat);

% error
error5 = e5_refit-e5_hat';
error6 = e6_refit-e6_hat';
disp('error5:'); disp(max(abs(error5(:,[3,6]))))
disp('error6:'); disp(max(abs(error6(:,[3,6]))))

%% 保存周期解
% save('../../subF_eom(CR3BP)/generalSolFFT_12', 'con', 'para', 'coe')
