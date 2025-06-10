set(0,'defaultAxesFontName', 'TimesSimSun');%������
set(0,'defaultTextFontName', 'TimesSimSun');%����

load('InitialValue_Polyfit.mat')

figure(1)
subplot(1,2,1)
plot(T_DRO_sample,J_sample,'o','MarkerSize',5)
hold on 
plot(T_DRO_all,J_all,'LineWidth',1.5)
hold off
legend('sample points')
xlabel('DRO������ת���ڱ�')
ylabel('�������')
set(gca,'FontSize',13)
    
subplot(1,2,2)
plot(T_DRO_sample,x0_sample,'o','MarkerSize',5)
hold on 
plot(T_DRO_all,x0_all,'LineWidth',1.5)
hold off
legend('sample points')
xlabel('DRO������ת���ڱ�')
ylabel('�����Ӽ�������ʱ��xֵ')
set(gca,'FontSize',13)