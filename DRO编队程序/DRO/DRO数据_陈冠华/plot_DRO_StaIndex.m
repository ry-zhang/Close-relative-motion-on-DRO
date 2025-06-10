%% 三体轨道中的线性化相对运动
% 2019-12-28
% by Yang Chihang
% email: ychhtl@foxmail.com
% close all
clear
addpath('../subF_eom(CR3BP)')
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
%% 常数与变量

load('data_DRO_')
T_norma_day = 4.34245177851179;
%% plot
figure(1)
v_StaInd = dataDRO(:,13:15);
v_StaInd = v_StaInd.*(abs(v_StaInd-1)>1e-5)+ones(size(v_StaInd)).*(abs(v_StaInd-1)<=1e-5);
v_StaInd(370:end,:) = 1;
plot(dataDRO(:,7)*T_norma_day,v_StaInd,'LineWidth',1.5,'Color',[0, 97, 182]/255); hold on
xlabel('周期 [day]'); ylabel('\itv');
set(gca,'FontSize',13)
title('DRO的稳定性指数'); 
grid on; grid minor
xlim([min(dataDRO(:,7)*T_norma_day), max(dataDRO(:,7)*T_norma_day)]) 
ylim(1+[-2e-4,4*1e-4])
% ylim([y_middle - x_ratio/2*y2x_ratio*x_diff, y_middle + x_ratio/2*y2x_ratio*x_diff]) 

% set(gcf,'Color',[255,255,255]/255);
% export_fig DRO.png -r600