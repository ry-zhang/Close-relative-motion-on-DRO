%% 载入数据
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
addpath('../../subF_eom(CR3BP)')
data_4_1 = load('FloquetEig_all4-1.mat');
data_4_2 = load('FloquetEig_all4-2.mat');

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);

%% 画图
T0_all = [data_4_1.J_period_all(:,2); data_4_2.J_period_all(:,2)]*con.T_norma_day*2*pi; % day
Meigva_all = [data_4_1.Meigva_all; data_4_2.Meigva_all];
T0_all = [T0_all(1:191);T0_all(191);T0_all(192:end)];
Meigva_all = [Meigva_all(1:191,:); [Meigva_all(191,1:4),1,1]; Meigva_all(192:end,:)];
theta = linspace(0,2*pi,1000);
x = cos(theta); y = sin(theta);

%------------------平面复特征根----------------------
figure(1)
subplot(2,1,1)
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
[~,index] = min(real(Meigva_all(:,5)));
patch([real(Meigva_all(index+1:end,5))'  NaN],[imag(Meigva_all(index+1:end,5))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(index+1:end,6))'  NaN],[imag(Meigva_all(index+1:end,6))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([3*real(Meigva_all(1:index,5))'  NaN],[3*imag(Meigva_all(1:index,5))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
xlabel('实部'); ylabel('虚部')
axis equal; set(gca,'FontSize',15)
xlim([-1.5,1.5]); ylim([-1.5,1.5]); 
% colorbar;  % Add a colorbar
hold off

subplot(2,1,2)
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
patch([real(Meigva_all(1:index,5))'  NaN],[imag(Meigva_all(1:index,5))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(1:index,6))'  NaN],[imag(Meigva_all(1:index,6))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([3*real(Meigva_all(index+1:end,5))'  NaN],[3*imag(Meigva_all(index+1:end,5))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
xlabel('实部'); ylabel('虚部')
axis equal; set(gca,'FontSize',15)
xlim([-1.5,1.5]); ylim([-1.5,1.5]); 
% colorbar;  % Add a colorbar
hold off

% 生成子图1
axes1 = axes('Position',[0.72,0.6467,0.16,0.27]); % 生成子图 
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
patch([real(Meigva_all(index+1:end,5))'  NaN],[imag(Meigva_all(index+1:end,5))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(index+1:end,6))'  NaN],[imag(Meigva_all(index+1:end,6))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
axis equal; set(gca,'FontSize',13)
xlim([0.9,1.1]); ylim([-0.1,0.1]); % 设置坐标轴范围

% 生成子图2
axes2 = axes('Position',[0.72,0.1730,0.16,0.27]); % 生成子图 
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
patch([real(Meigva_all(1:index,5))'  NaN],[imag(Meigva_all(1:index,5))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(1:index,6))'  NaN],[imag(Meigva_all(1:index,6))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
axis equal; set(gca,'FontSize',13)
xlim([0.9,1.1]); ylim([-0.1,0.1]); % 设置坐标轴范围

exportgraphics(gcf,'NormalEigenVa.jpg','Resolution',600)

%------------------法向复特征根----------------------
figure(2)
subplot(2,1,1)
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
[~,index] = min(real(Meigva_all(:,1)));
patch([real(Meigva_all(index+1:end,1))'  NaN],[imag(Meigva_all(index+1:end,1))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(index+1:end,2))'  NaN],[imag(Meigva_all(index+1:end,2))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([3*real(Meigva_all(1:index,2))'  NaN],[3*imag(Meigva_all(1:index,2))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
xlabel('实部'); ylabel('虚部')
axis equal; set(gca,'FontSize',15)
xlim([-1.5,1.5]); ylim([-1.5,1.5]); 
% colorbar;  % Add a colorbar
hold off

subplot(2,1,2)
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
patch([real(Meigva_all(1:index,1))'  NaN],[imag(Meigva_all(1:index,1))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(1:index,2))'  NaN],[imag(Meigva_all(1:index,2))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([3*real(Meigva_all(index+1:end,2))'  NaN],[3*imag(Meigva_all(index+1:end,2))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
xlabel('实部'); ylabel('虚部')
axis equal; set(gca,'FontSize',15)
xlim([-1.5,1.5]); ylim([-1.5,1.5]); 
% colorbar;  % Add a colorbar
hold off

% 生成子图
axes3 = axes('Position',[0.72,0.6467,0.16,0.27]); % 生成子图 
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
patch([real(Meigva_all(index+1:end,1))'  NaN],[imag(Meigva_all(index+1:end,1))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(index+1:end,2))'  NaN],[imag(Meigva_all(index+1:end,2))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
axis equal; set(gca,'FontSize',13)
xlim([0.9,1.1]); ylim([-0.1,0.1]); % 设置坐标轴范围
exportgraphics(gcf,'PlanarEigenVa.jpg','Resolution',600)

figure(3)
plot(x+y*1i,'LineWidth',1,'Color',0.7*[1,1,1]); hold on;
patch([real(Meigva_all(1:index,1))'  NaN],[imag(Meigva_all(1:index,1))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([real(Meigva_all(1:index,2))'  NaN],[imag(Meigva_all(1:index,2))' NaN],[T0_all(1:index)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')
patch([3*real(Meigva_all(index+1:end,2))'  NaN],[3*imag(Meigva_all(index+1:end,2))' NaN],[T0_all(index+1:end)' NaN],'LineWidth',3,'EdgeColor','interp','FaceColor','none')

xlabel('实部'); ylabel('虚部')
axis equal; set(gca,'FontSize',15)
xlim([-1.5,1.5]); ylim([-1.5,1.5]); 
colorbar;  % Add a colorbar
ylabel(colorbar,'周期 [day]','FontSize',15);
hold off
exportgraphics(gcf,'Colorbar.jpg','Resolution',600)