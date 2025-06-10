load('NaturalSSF_DRO1_90day_DE.mat')
plot(MaxDist_all,'LineWidth',1.5); hold on
load('NaturalSSF_DRO1_90day_GA.mat')
plot(MaxDist_all,'LineWidth',1.5); hold off
set(gca,'FontSize',13)
xlabel('优化次数')
ylabel('优化结果')
% ylabel('优化结果(距参考轨道的最大距离)[km]')
legend('朱小龙 DE','MATLBA GA')