%% CW方程的xy方向振动奇异点分析
% 2022-3-9
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear

format longg
format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultAxesFontSize', 13);%坐标轴
set(0,'defaultTextFontSize', 13);%文字
set(0,'defaultLineLineWidth',1.5)

%% 计算t0与dt对应的phi21的值
length_dt = 1001;
dt_sample = linspace(0,3,length_dt);
phidet_all = zeros(length_dt,1);
n = 2*pi;

for jj_loop = 1:length_dt
    dt = dt_sample(jj_loop);
    Phi_tf = PhiCWxy(dt,n);
%     phidet_all(jj_loop) = det(Phi_tf([1,2],[3,4]));
    phidet_all(jj_loop) = (8-8*cos(n*dt)-3*n*dt*sin(n*dt))/n^2;
end

%% plot
% ---------------------------t0,dt与||\Phi_{xy}^{v}||----------------------------------
figure(1)
plot(dt_sample,phidet_all')

xlabel('\Delta{\itt} [{\itT}_0]'); ylabel('$$||\Phi_{xy}^{v}||$$','Interpreter','latex')
grid on; grid minor; box on;
