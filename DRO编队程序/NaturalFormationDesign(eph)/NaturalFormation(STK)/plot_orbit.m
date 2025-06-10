set(0,'defaultAxesFontName', 'TimesSimSun','defaultTextFontName', 'TimesSimSun');
% 关于TimesSimSun字体文件请参看博文https://blog.csdn.net/qq_15950515/article/details/122991780
set(0,'defaultAxesFontSize',15,'defaultTextFontSize',15)
set(0,'defaultLineLineWidth',1.5)

load('DRO_SSF.mat')
transferSFF = auxSFF.transferSFF;
naturalSFF = auxSFF.naturalSFF;
sepSFF = auxSFF.sepSFF;

ii_forma = 1;

%% ---------月心地月旋转系LVLH-----------
xx_MCEMRLVLH_REL_sep = sepSFF.xx_MCEMRLVLH_REL;
xx_MCEMRLVLH_REL_trans = transferSFF(ii_forma).xx_MCEMRLVLH_REL;
xx_MCEMRLVLH_REL_natural =  naturalSFF(ii_forma).xx_MCEMRLVLH_REL;

% 画图
figure(1)
p0 = plot3(0,0,0,'ks','MarkerSize',5); hold on;
p1 = plot3(xx_MCEMRLVLH_REL_trans(:,1), xx_MCEMRLVLH_REL_trans(:,2), xx_MCEMRLVLH_REL_trans(:,3),'r','LineWidth',1.5); 
p2 = plot3(xx_MCEMRLVLH_REL_natural(:,1), xx_MCEMRLVLH_REL_natural(:,2), xx_MCEMRLVLH_REL_natural(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5);
p3 = plot3(xx_MCEMRLVLH_REL_natural(1,1), xx_MCEMRLVLH_REL_natural(1,2), xx_MCEMRLVLH_REL_natural(1,3),'g^');
p4 = plot3(xx_MCEMRLVLH_REL_natural(end,1), xx_MCEMRLVLH_REL_natural(end,2), xx_MCEMRLVLH_REL_natural(end,3),'rv');
if ii_forma == 1
    p01 = plot3(xx_MCEMRLVLH_REL_sep(:,1), xx_MCEMRLVLH_REL_sep(:,2), xx_MCEMRLVLH_REL_sep(:,3),'Color',[16, 132, 254]/255,'LineWidth',1.5); 
    legend([p0,p01,p1,p2,p3,p4],{'主星','分离轨道','转移轨道','任务轨道','初始位置','终端位置'},'Location','northeastoutside');
    xx_MCEMRLVLH_REL_all = [xx_MCEMRLVLH_REL_sep; xx_MCEMRLVLH_REL_trans; xx_MCEMRLVLH_REL_natural];
else
    legend([p0,p1,p2,p3,p4],{'主星','转移轨道','任务轨道','初始位置','终端位置'},'Location','northeastoutside');
    xx_MCEMRLVLH_REL_all = [xx_MCEMRLVLH_REL_trans; xx_MCEMRLVLH_REL_natural];
end

box on; grid on; grid minor; 
% hold off; 
axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); zlabel('\itz_L \rm[km]')
xlim(1.5*[min(xx_MCEMRLVLH_REL_all(:,1)),max(xx_MCEMRLVLH_REL_all(:,1))]); 
ylim(1.2*[min(xx_MCEMRLVLH_REL_all(:,2)),max(xx_MCEMRLVLH_REL_all(:,2))]); 
set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
title('月心地月旋转系LVLH')
view(0,90)
set(gcf,'Color',[255,255,255]/255);
% print(gcf, '-dpng', 'MCRLVLH.png','-r600');

%% ---------地心J2000 LVLH-----------
xx_j2kLVLH_REL_sep = sepSFF.xx_j2kLVLH_REL;
xx_j2kLVLH_REL_trans = transferSFF(ii_forma).xx_j2kLVLH_REL;
xx_j2kLVLH_REL_natural =  naturalSFF(ii_forma).xx_j2kLVLH_REL;

% 画图
figure(2)
p0 = plot3(0,0,0,'ks','MarkerSize',5); hold on
p1 = plot3(xx_j2kLVLH_REL_trans(:,1),xx_j2kLVLH_REL_trans(:,2),xx_j2kLVLH_REL_trans(:,3),'r','LineWidth',1.5);
p2 = plot3(xx_j2kLVLH_REL_natural(:,1), xx_j2kLVLH_REL_natural(:,2), xx_j2kLVLH_REL_natural(:,3),'Color',[237, 177, 32]/255,'LineWidth',1.5);
p3 = plot3(xx_j2kLVLH_REL_natural(1,1), xx_j2kLVLH_REL_natural(1,2), xx_j2kLVLH_REL_natural(1,3),'g^');
p4 = plot3(xx_j2kLVLH_REL_natural(end,1), xx_j2kLVLH_REL_natural(end,2), xx_j2kLVLH_REL_natural(end,3),'rv');
if ii_forma == 1
    p01 = plot3(xx_j2kLVLH_REL_sep(:,1), xx_j2kLVLH_REL_sep(:,2), xx_j2kLVLH_REL_sep(:,3),'Color',[16, 132, 254]/255,'LineWidth',1.5); 
    legend([p0,p01,p1,p2,p3,p4],{'主星','分离轨道','转移轨道','任务轨道','初始位置','终端位置'},'Location','northeastoutside');
    xx_j2kLVLH_REL_all = [xx_j2kLVLH_REL_sep; xx_j2kLVLH_REL_trans; xx_j2kLVLH_REL_natural];
else
    legend([p0,p1,p2,p3,p4],{'主星','转移轨道','任务轨道','初始位置','终端位置'},'Location','northeastoutside');
    xx_j2kLVLH_REL_all = [xx_j2kLVLH_REL_trans; xx_j2kLVLH_REL_natural];
end
% hold off; 
box on; grid on; grid minor; 
axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); zlabel('\itz_L \rm[km]')
xlim(1.2*[min(xx_j2kLVLH_REL_all(:,1)),max(xx_j2kLVLH_REL_all(:,1))]); 
ylim(1.1*[min(xx_j2kLVLH_REL_all(:,2)),max(xx_j2kLVLH_REL_all(:,2))]); 
title('地心J2000LVLH')
set(gca,'FontSize',15);
view(0,90)

set(gcf,'Color',[255,255,255]/255);
print(gcf, '-dpng', 'ECILVLH.png','-r600');
