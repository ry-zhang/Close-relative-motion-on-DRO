set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

load('InitialValue_Polyfit.mat')

figure(1)
subplot(1,2,1)
plot(T_DRO_sample,J_sample,'o','MarkerSize',5)
hold on 
plot(T_DRO_all,J_all,'LineWidth',1.5)
hold off
legend('sample points')
xlabel('DRO与月球公转周期比')
ylabel('轨道能量')
set(gca,'FontSize',13)
    
subplot(1,2,2)
plot(T_DRO_sample,x0_sample,'o','MarkerSize',5)
hold on 
plot(T_DRO_all,x0_all,'LineWidth',1.5)
hold off
legend('sample points')
xlabel('DRO与月球公转周期比')
ylabel('穿过庞加莱截面时的x值')
set(gca,'FontSize',13)