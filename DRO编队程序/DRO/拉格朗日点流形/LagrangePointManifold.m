
% 拉格朗日点的流形结构
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
% 参考文献：刘林，胡松杰，王歆，航天器动力学引论,章节7.4
clear
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultLineLineWidth',2)

format longg
format compact
opts = odeset('RelTol',1e-13,'AbsTol',1e-20);

mu = 0.01215;
L = lagrange_points(mu);
t_sample = linspace(0,2*pi,1000);
%% 共线拉格朗日点
% 线性化相对运动方程的系数
% 以L2为例
r1 = abs(L(2,1)-(-mu)); 
r2 = abs(L(2,1)-(1-mu));
C0 = (1-mu)/r1^3 + mu/r2^3;
Omegaxx = 1+2*C0; 
Omegayy = 1-C0;
Omegazz = -C0;
Omegaxy = 0;

A21 = [Omegaxx, Omegaxy, 0;
    Omegaxy, Omegayy, 0;
    0,0,Omegazz];
A22 = [0,2,0; -2,0,0; 0,0,0];
A = [zeros(3), eye(3);
    A21, A22];
[eigve,eigva] = eig(A);
eigva = eig(A);
% sol3 = ode113(@(t,x)lagrangeRelMot(t,x,A),[0 dt], [para(i_index).x0_DRO, EigenVector_all(i_index).p3'], opts);
lambda1 = eigva(1); lambda3 = eigva(3); lambda5 = eigva(5);
d1 = lambda1; d2 = imag(lambda3); d3 = imag(lambda5);
k1 = (d1^2-1-2*C0)/2/d1; k2 = (d2^2+1+2*C0)/2/d2;

r_unstable = 1e-5*[exp(d1*t_sample); k1*exp(d1*t_sample); zeros(size(t_sample))];
r_stable = 7*[exp(-d1*t_sample); -k1*exp(-d1*t_sample); zeros(size(t_sample))];
r_period1 = 0.3*[cos(d2*t_sample); -k2*sin(d2*t_sample); zeros(size(t_sample))];
r_period2 = 0.3*[sin(d2*t_sample); +k2*cos(d2*t_sample); zeros(size(t_sample))];
r_period_z1 = 0.7*[zeros(size(t_sample)); zeros(size(t_sample)); cos(d3*t_sample)];
r_period_z2 = 0.7*[zeros(size(t_sample)); zeros(size(t_sample)); sin(d3*t_sample)];

% 画图
figure(1);
color_all = get(gca,'colororder');
p1 = plot3(r_stable(1,:),r_stable(2,:),r_stable(3,:),'Color',color_all(1,:)); hold on;
p1 = plot3(-r_stable(1,:),-r_stable(2,:),-r_stable(3,:),'Color',color_all(1,:));
p2 = plot3(r_unstable(1,:),r_unstable(2,:),r_unstable(3,:),'Color',color_all(2,:)); 
p2 = plot3(-r_unstable(1,:),-r_unstable(2,:),-r_unstable(3,:),'Color',color_all(2,:)); 
p3 = plot3(r_period1(1,:),r_period1(2,:),r_period1(3,:),'Color',color_all(3,:));
p3 = plot3(r_period2(1,:),r_period2(2,:),r_period2(3,:),'Color',color_all(3,:));
p4 = plot3(r_period_z1(1,:),r_period_z1(2,:),r_period_z1(3,:),'Color',color_all(4,:));
p4 = plot3(r_period_z2(1,:),r_period_z2(2,:),r_period_z2(3,:),'Color',color_all(4,:));
pgon = polyshape(3*[-1,1; 1,1; 1,-1; -1,-1]);
p5 = plot(pgon,'FaceColor',0.5*[1,1,1],'FaceAlpha',0.5);
hold off
legend([p5,p1,p2,p3,p4],{'地月绕转平面','面内稳定流形','面内不稳定流形','面内中心流形','法向中心流形'},...
    'Location','northeast')
grid on; box on
xlim([-2,2]); ylim([-1.2,1.2]); zlim(0.7*[-1,1])
set(gca,'FontSize',15)
xlabel('\itx'); ylabel('\ity'); zlabel('\itz');

view([-15,30])
exportgraphics(gcf,'L2Manifold.jpg','Resolution',600)

%% 三角拉格朗日点
% 线性化相对运动方程的系数
% 以L2为例
r1 = abs(L(4,1)-(-mu)); 
r2 = abs(L(4,1)-(1-mu));
C0 = (1-mu)/r1^3 + mu/r2^3;
Omegaxx = 3/4; 
Omegayy = 9/4;
Omegazz = -C0;
Omegaxy = 3*sqrt(3)/2*(mu-1/2);

A21 = [Omegaxx, Omegaxy, 0;
    Omegaxy, Omegayy, 0;
    0,0,Omegazz];
A22 = [0,2,0; -2,0,0; 0,0,0];
A = [zeros(3), eye(3);
    A21, A22];
[eigve,eigva] = eig(A);
% eigva = eig(A);
% eigva = sort(eigva);
[~,r_center1] = ode113(@(t,x)lagrangeRelMot(t,x,A),[0 10*pi], 1.5*real(eigve(:,1)), opts); r_center1 = r_center1';
[~,r_center2] = ode113(@(t,x)lagrangeRelMot(t,x,A),[0 10*pi], 1.5*real(eigve(:,3)), opts); r_center2 = r_center2';
[~,r_center_z] = ode113(@(t,x)lagrangeRelMot(t,x,A),[0 10*pi], 3*real(eigve(:,5)), opts); r_center_z = r_center_z';
% lambda1 = eigva(1); lambda3 = eigva(3); lambda5 = eigva(5);
% d1 = lambda1; d2 = imag(lambda3); d3 = imag(lambda5);
% k1 = (d1^2-1-2*C0)/2/d1; k2 = (d2^2+1+2*C0)/2/d2;
% 
% r_unstable = 1e-5*[exp(lambda1*t_sample); k1*exp(lambda1*t_sample); zeros(size(t_sample))];
% r_stable = 7*[exp(-lambda1*t_sample); -k1*exp(-lambda1*t_sample); zeros(size(t_sample))];
% r_period1 = 0.3*[cos(d2*t_sample); -k2*sin(d2*t_sample); zeros(size(t_sample))];
% r_period2 = 0.3*[sin(d2*t_sample); +k2*cos(d2*t_sample); zeros(size(t_sample))];
% r_period_z1 = 0.7*[zeros(size(t_sample)); zeros(size(t_sample)); cos(d3*t_sample)];
% r_period_z2 = 0.7*[zeros(size(t_sample)); zeros(size(t_sample)); sin(d3*t_sample)];

% 画图
figure(2);
color_all = get(gca,'colororder');
p1 = plot3(r_center1(1,:),r_center1(2,:),r_center1(3,:),'Color',color_all(3,:)); hold on;
p1 = plot3(-r_center1(1,:),-r_center1(2,:),-r_center1(3,:),'Color',color_all(3,:));
p2 = plot3(r_center2(1,:),r_center2(2,:),r_center2(3,:),'Color',color_all(5,:)); 
p2 = plot3(-r_center2(1,:),-r_center2(2,:),-r_center2(3,:),'Color',color_all(5,:)); 
p3 = plot3(r_center_z(1,:),r_center_z(2,:),r_center_z(3,:),'Color',color_all(4,:));
p3 = plot3(-r_center_z(1,:),-r_center_z(2,:),-r_center_z(3,:),'Color',color_all(4,:));
pgon = polyshape(3*[-1,1; 1,1; 1,-1; -1,-1]);
p4 = plot(pgon,'FaceColor',0.5*[1,1,1],'FaceAlpha',0.5);
hold off
legend([p4,p1,p2,p3],{'地月绕转平面','面内中心流形(短周期)','面内中心流形(长周期)','法向中心流形'},...
    'Location','northeast')
grid on; box on
xlim([-2,2]); ylim([-1.2,1.2]); zlim(0.7*[-1,1])
set(gca,'FontSize',15)
xlabel('\itx'); ylabel('\ity'); zlabel('\itz');

view([-15,30])
exportgraphics(gcf,'L5Manifold.jpg','Resolution',600)