clear
close all
format longg
format compact
warning on
dbstop if error

set(0,'defaultAxesFontName', 'TimesSimSun','defaultTextFontName', 'TimesSimSun');
set(0,'defaultAxesFontSize',13,'defaultTextFontSize',13)
set(0,'defaultLineLineWidth',1.5)

%% 
iDRO = 0;
scale_all = [1,10,100,1000]; % km
for ii_index = 1:4
    scale = scale_all(ii_index);
    load(['2DFormation_DR0',num2str(iDRO),'_Scale',num2str(scale)]);
    f1 = figure(312); set(gca,'YScale','log');
%     dv_eph_all = dv_eph_all.*(r_LVLH_rel_eph_mane_error_all<1e-3) + dv_CRTBP_all*ones(1,size(dv_eph_all,2)).*(r_LVLH_rel_eph_mane_error_all>=1e-3);
    
    r_error_all_filt = max(r_LVLH_rel_eph_mane_error_all(:,2:2:20), r_LVLH_rel_eph_mane_error_all(:,3:2:21))>=1e-3;
    dv_eph_all_filt = dv_eph_all(:,2:2:20) + dv_eph_all(:,3:2:21);
    for jj_index = 1:size(dv_eph_all_filt,1)
        dv_eph_all_filt(jj_index,:) = dv_eph_all_filt(jj_index,:).*(~r_error_all_filt(jj_index,:)) + ...
            mean(dv_eph_all_filt(jj_index,~r_error_all_filt(jj_index,:))).*r_error_all_filt(jj_index,:);
    end
%     boxplot(dt_all/para.T0,dv_eph_all');
%     boxplot(dv_eph_all_filt(1:end,:)'); hold on
    
%     plot(dt_all/para.T0,max(dv_eph_all_filt,[],2),'.','Color',[0,0.4470,0.7410]*1.2);
%     plot(dt_all/para.T0,min(dv_eph_all_filt,[],2),'.','Color',[0,0.4470,0.7410]*1.2);
    p_pat = patch([dt_all,dt_all(end:-1:1)]/para.T0,[max(dv_eph_all_filt,[],2)',min(dv_eph_all_filt(end:-1:1,:),[],2)'],[0.938, 0.814, 0.565]);
    p_pat.EdgeAlpha = 0;
    hold on
    p0 = plot(dt_all/para.T0,2*dv_CRTBP_all,'Color','k');
    p1 = plot(dt_all/para.T0,mean(dv_eph_all_filt,2),'Color',[0.925, 0.625, 0.025]);  % [237, 177, 32]/255
end
h1 = figure(312);
ylim([1e-3,2e3])
xlabel('\Delta{\itt} [{\itT}_0]'); ylabel('\Delta{\itv} [m/s]')
legend([p0,p1,p_pat],{'CRTBP','星历平均值','星历上下界'},'Location','north')
% legend([p0,p1,p_pat],{'CRTBP','eph average','eph bounds'},'Location','north')
% xlabel('\Delta{\itt} [{\itT}]');
grid on; box on
% 标上不同尺度的轨道

% 搞个放大图
iaxes = MagInset(h1, -1, [0.28 0.34 1.8 3.3], [0.58 0.77 10 270], {'NE','NW';'SE','SW'});
set(gca,'fontsize',13);
% iaxes = MagInset(h1, -1, [0.25 0.35 1.5,2.15], [0.58 0.77 4 200], {'NE','NW';'SE','SW'});
hold off

set(gcf,'Color',[255,255,255]/255);
export_fig 2DForma_dv.png -r600
