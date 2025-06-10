%% 高精度模型下相对运动变轨计算程序
% v0：杨驰航(2022/02/22)，email: ychhtl@foxmail.com
%     创建并撰写核心程序

%% 程序说明
% [Dv_j2kLVLH, dv12, rf_error, minDist] = 
% NaturalSFFTransfer(auxSFF,tspan_transfer,ii_segment,dv0,isplot,DynamicModel)
% 
% Input arguments:
% -------------------------------------------------------------------------
% auxSFF.jd0:                           DRO任务入轨位置(TBD儒略日)
%       .x00_j2k_target:                DRO任务入轨位置轨道状态
%       .naturalSFF:                    储存三次编队轨道数据的结构体
%       .naturalSFF(i).scale_ref:       第i个编队的轨道尺度
%       .naturalSFF(i).t0:              第i个编队的优化起始位置(相对于jd0的秒数)
%       .naturalSFF(i).tf:              第i个编队的优化终端位置(相对于jd0的秒数)
%       .naturalSFF(i).x0_MCRLVLH_rel:  t0处的MCRLVLH下的相对运动状态
%       .naturalSFF(i).x0_MCR_target:   t0处的MCR下的主星状态
%       .naturalSFF(i).x0_j2kLVLH_rel:  t0处的j2kLVLH下的相对运动状态
%       .naturalSFF(i).x0_j2k_target:   t0处的j2k下的主星状态
% tspan_transfer:                       转移始末时间[tm0,tmf]
% ii_segment:                           优化转移轨道段序号(1:分离-段1, 2:段1-段2, 3:段2-段3)
% dv0:                                  分离脉冲，仅在ii_segment为1时有效
% isplot:                               画图标志位(0/1)
% DynamicModel:                         动力学模型('SEMephemeris'或'HighFidelity')
% 
% Output arguments:
% -------------------------------------------------------------------------
% Dv12_j2kLVLH:      naturalSFF(i).t0/tf 处的变轨处脉冲 [2*3] (m/s,第一行Dv1,第二行Dv2)
% dv12_j2kLVLH:      naturalSFF(i).t0/tf 处的变轨处脉冲大小 [1*2] (m/s,1为||Dv1||,2为||Dv2||)
% rf_error:           naturalSFF(i).tf 处的末端位置误差 [1*1] (km)
%                    如果rf_error大于1m，则Dv12_j2kLVLH与dv12_j2kLVLH的返回值是无穷
% minDist:           naturalSFF(i).t0-tf 转移轨道上的主副星最小距离 [1*1] (km)
% 
% 程序中航天器状态变量命名规则
% -------------------------------------------------------------------------
% 命名格式: 变量名_坐标系_航天器_后缀
% 变量名:   x0-初始状态，xf-末端状态，xx-整个轨道状态，对于a/r/v，采用同样命名规则
% 坐标系:   MCR-月心旋转系，j2k-地心J2000坐标系，MCRLVLH-MCR下的LVLH相对系，j2kLVLH-j2k下的LVLH相对系
% 航天器:   target-DRO目标航天器，chaser-绕飞航天器，rel-两者相对运动
% 后缀:     可用于标识额外信息，可有多段。例如描述不同阶段：Last\Next\transfer,于描述不同模型：CRTBP/eph
% 例:       xx_MCR_chaser_Next \ x0_MCRLVLH_rel_Last \ xf_j2kLVLH_rel_Last
% 
%% 主程序
function auxSFF_dv = NaturalSFFTransfer(auxSFF,tspan_transfer,ii_segment,dv0,isplot,DynamicModel)

naturalSFF = auxSFF.naturalSFF;
x00_j2k_target = auxSFF.x00_j2k_target;
tm0 = tspan_transfer(1);
tmf = tspan_transfer(2);
dt_transfer = diff(tspan_transfer);

%% 星历初始化
aux = []; 
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.jd0 = auxSFF.jd0;

if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,1); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux = initialize(aux); % 初始化
else
    error('Wrong DynamicModel')
end
% x00_MCR_target = T_ECJ2k2Rot(aux.jd0,x00_j2k_target,[],aux.C_Mat, 'MCEMR',0);

%% 变轨优化
if ii_segment == 1
    % Last segment（分离滑行阶段）
    epoch_SSF_Last = naturalSFF(ii_segment).t0;
    tspan_sec_Last = [epoch_SSF_Last, tm0];
    t_sample_Last = linspace(tspan_sec_Last(1),tspan_sec_Last(2),max(500,ceil(diff(tspan_sec_Last)/3600/4)+1));
    % 主星: 由分离位置(epoch_SSF_Last)递推至分离结束位置(也即变轨开始位置tm0)
    x0_MCR_target_Last = naturalSFF(ii_segment).x0_MCR_target;
    [xx_MCR_target_Last,aa_MCR_target_Last] = Propagate_EphRotFrame(x0_MCR_target_Last,tspan_sec_Last,t_sample_Last,aux);
    % 副星: 由Last segment初始位置递推至Last segment末端位置(也即变轨初始位置)
    dv0_j2kLVLH = 1e-3*dv0; % m/s -> km/s
    xx_j2k_target_Last = T_Rot2ECJ2k(aux.jd0+t_sample_Last/86400,xx_MCR_target_Last,aux.C_Mat, 'MCEMR',0);
    aa_j2k_target_Last = eomj2kMtx(xx_j2k_target_Last,t_sample_Last,aux);
    x0_j2kLVLH_rel_Last = [0,0,0,dv0_j2kLVLH];
    x0_j2k_rel_Last = T_TCO2TCR_eph(x0_j2kLVLH_rel_Last,xx_j2k_target_Last(1,:),aa_j2k_target_Last(1,:),'LVLH');
    x0_j2k_chaser_Last = xx_j2k_target_Last(1,:) + x0_j2k_rel_Last;
    x0_MCR_chaser_Last = T_ECJ2k2Rot(aux.jd0+epoch_SSF_Last/86400,x0_j2k_chaser_Last,[],aux.C_Mat, 'MCEMR',0);

    xx_MCR_chaser_Last = Propagate_EphRotFrame(x0_MCR_chaser_Last,tspan_sec_Last,t_sample_Last,aux);
    xx_MCR_rel_Last = xx_MCR_chaser_Last-xx_MCR_target_Last;
    xx_MCRLVLH_rel_Last = T_TCR2TCO_eph(xx_MCR_rel_Last,xx_MCR_target_Last,aa_MCR_target_Last,'LVLH');
    xf_MCRLVLH_rel_Last = xx_MCRLVLH_rel_Last(end,:);
    xf_MCR_target_Last = xx_MCR_target_Last(end,:);
else
    % Last segment
    epoch_SSF_Last = naturalSFF(ii_segment-1).t0;
    tspan_sec_Last = [epoch_SSF_Last, tm0];
    t_sample_Last = linspace(tspan_sec_Last(1),tspan_sec_Last(2),max(1500,ceil(diff(tspan_sec_Last)/3600/4)+1));
    % 主星: 由Last segment初始位置递推至Last segment末端位置(也即变轨初始位置)
    x0_MCR_target_Last = naturalSFF(ii_segment-1).x0_MCR_target;
    [xx_MCR_target_Last,aa_MCR_target_Last] = Propagate_EphRotFrame(x0_MCR_target_Last,tspan_sec_Last,t_sample_Last,aux);
    % 副星: 由Last segment初始位置递推至Last segment末端位置(也即变轨初始位置)
    x0_MCRLVLH_rel_Last = naturalSFF(ii_segment-1).x0_MCRLVLH_rel;
    x0_MCR_chaser_Last = T_TCO2TCR_eph(x0_MCRLVLH_rel_Last,x0_MCR_target_Last,aa_MCR_target_Last(1,:),'LVLH')+x0_MCR_target_Last;
    xx_MCR_chaser_Last = Propagate_EphRotFrame(x0_MCR_chaser_Last,tspan_sec_Last,t_sample_Last,aux);
    xx_MCR_rel_Last = xx_MCR_chaser_Last-xx_MCR_target_Last;
    xx_MCRLVLH_rel_Last = T_TCR2TCO_eph(xx_MCR_rel_Last,xx_MCR_target_Last,aa_MCR_target_Last,'LVLH');
    xf_MCRLVLH_rel_Last = xx_MCRLVLH_rel_Last(end,:);
    xf_MCR_target_Last = xx_MCR_target_Last(end,:);
end

% Next segment
epoch_SSF_Next = naturalSFF(ii_segment).t0;
tspan_sec_Next = [epoch_SSF_Next, tmf];
t_sample_Next = linspace(tspan_sec_Next(1),tspan_sec_Next(2),max(500,ceil(diff(tspan_sec_Next)/3600/4)+1));
% 主星: 由Next segment初始位置(也即变轨初始位置)递推至变轨末端位置
x0_MCR_target_Next = naturalSFF(ii_segment).x0_MCR_target;
[xx_MCR_target_Next,af_MCR_target_Next] = Propagate_EphRotFrame(x0_MCR_target_Next,tspan_sec_Next,t_sample_Next,aux);
% 副星: 由Next segment初始位置(也即变轨初始位置)递推至变轨末端位置
x0_MCRLVLH_rel_Next = naturalSFF(ii_segment).x0_MCRLVLH_rel;
x0_MCR_chaser_Next = T_TCO2TCR_eph(x0_MCRLVLH_rel_Next,x0_MCR_target_Next,af_MCR_target_Next(1,:),'LVLH')+x0_MCR_target_Next;
xx_MCR_chaser_Next = Propagate_EphRotFrame(x0_MCR_chaser_Next,tspan_sec_Next,t_sample_Next,aux);
xx_MCR_rel_Next = xx_MCR_chaser_Next-xx_MCR_target_Next;
xx_MCRLVLH_rel_Next = T_TCR2TCO_eph(xx_MCR_rel_Next,xx_MCR_target_Next,af_MCR_target_Next,'LVLH');
xf_MCRLVLH_rel_Next = xx_MCRLVLH_rel_Next(end,:);

% 变轨迭代
rv0f_MCRLVLH_rel = [xf_MCRLVLH_rel_Last; xf_MCRLVLH_rel_Next]; % 变轨起始与终端位置
aux_transfer = aux;
aux_transfer.jd0 = aux.jd0+tm0/86400;
aux_transfer.t0UTC = [];
if strcmp(DynamicModel,'HighFidelity')
    aux_transfer = SPICEinitialize(aux_transfer,1); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux_transfer = initialize(aux_transfer); % 初始化
end
[~,rv0f_MCRLVLH_rel_mane,~,~,~,rf_error] = ...
    forcedRelMotion_eph(xf_MCR_target_Last,rv0f_MCRLVLH_rel,dt_transfer,aux_transfer,'LVLH',0,0);
Dv1_MCRLVLH = (rv0f_MCRLVLH_rel_mane(1,4:6) - rv0f_MCRLVLH_rel(1,4:6))*1e3;
Dv2_MCRLVLH = (rv0f_MCRLVLH_rel(2,4:6) - rv0f_MCRLVLH_rel_mane(2,4:6))*1e3;
Dv_MCRLVLH = [Dv1_MCRLVLH;Dv2_MCRLVLH];

%% 根据优化结果计算LVLH坐标系下的相对运动
% Last segment轨道转化到j2k下
xx_j2k_chaser_Last = T_Rot2ECJ2k(aux.jd0+t_sample_Last/86400,xx_MCR_chaser_Last,aux.C_Mat, 'MCEMR',0);
xx_j2k_target_Last = T_Rot2ECJ2k(aux.jd0+t_sample_Last/86400,xx_MCR_target_Last,aux.C_Mat, 'MCEMR',0);
xx_j2k_rel_Last = xx_j2k_chaser_Last-xx_j2k_target_Last;
xx_j2kLVLH_rel_Last = T_TCR2TCO_eph(xx_j2k_rel_Last,xx_j2k_target_Last,[],'LVLH');

% Next segment轨道,变轨后递推
tspan_sec_Next2 = [tmf,naturalSFF(ii_segment).tf];
t_sample_Next2 = linspace(tspan_sec_Next2(1),tspan_sec_Next2(2),max(1500,ceil(diff(tspan_sec_Next2)/3600/4)+1));
% 主/副星: 由Next segment变轨末端位置递推至Next segment末端位置
[xx_MCR_target_Next2,af_MCR_target_Next2] = Propagate_EphRotFrame(xx_MCR_target_Next(end,:),tspan_sec_Next2,t_sample_Next2,aux);
[xx_MCR_chaser_Next2,~] = Propagate_EphRotFrame(xx_MCR_chaser_Next(end,:),tspan_sec_Next2,t_sample_Next2,aux);
xx_MCR_rel_Next2 = xx_MCR_chaser_Next2-xx_MCR_target_Next2;
xx_MCRLVLH_rel_Next2 = T_TCR2TCO_eph(xx_MCR_rel_Next2,xx_MCR_target_Next2,af_MCR_target_Next2,'LVLH');
% 主/副/相对运动: 转换至j2k
xx_j2k_target_Next2 = T_Rot2ECJ2k(aux.jd0+t_sample_Next2/86400,xx_MCR_target_Next2,aux.C_Mat, 'MCEMR',0);
xx_j2k_chaser_Next2 = T_Rot2ECJ2k(aux.jd0+t_sample_Next2/86400,xx_MCR_chaser_Next2,aux.C_Mat, 'MCEMR',0);
xx_j2k_rel_Next2 = xx_j2k_chaser_Next2-xx_j2k_target_Next2;
xx_j2kLVLH_rel_Next2 = T_TCR2TCO_eph(xx_j2k_rel_Next2,xx_j2k_target_Next2,[],'LVLH');

% Transfer轨道,变轨递推
tspan_sec_transfer = [tm0,tmf];
t_sample_transfer = tm0+linspace(0,dt_transfer,500);
% 主/副/相对运动: 由变轨初始位置递推至变轨末端位置,
[xx_MCR_target_transfer,aa_MCR_target_transfer] = Propagate_EphRotFrame(xf_MCR_target_Last,tspan_sec_transfer,t_sample_transfer,aux);
xf_MCR_chaser_Last = T_TCO2TCR_eph(rv0f_MCRLVLH_rel_mane(1,:),xf_MCR_target_Last,aa_MCR_target_transfer(1,:),'LVLH') + xf_MCR_target_Last;
xx_MCR_chaser_transfer = Propagate_EphRotFrame(xf_MCR_chaser_Last,tspan_sec_transfer,t_sample_transfer,aux);
xx_MCR_rel_transfer = xx_MCR_chaser_transfer-xx_MCR_target_transfer;
xx_MCRLVLH_rel_transfer = T_TCR2TCO_eph(xx_MCR_rel_transfer,xx_MCR_target_transfer,aa_MCR_target_transfer,'LVLH');

% 主/副/相对运动: 转换到j2k
xx_j2k_target_transfer = T_Rot2ECJ2k(aux.jd0+t_sample_transfer/86400,xx_MCR_target_transfer,aux.C_Mat, 'MCEMR',0);
xx_j2k_chaser_transfer = T_Rot2ECJ2k(aux.jd0+t_sample_transfer/86400,xx_MCR_chaser_transfer,aux.C_Mat, 'MCEMR',0);
xx_j2k_rel_transfer = xx_j2k_chaser_transfer-xx_j2k_target_transfer;
xx_j2kLVLH_rel_transfer = T_TCR2TCO_eph(xx_j2k_rel_transfer,xx_j2k_target_transfer,[],'LVLH');

% 将变轨脉冲转化到j2kLVLH坐标系下
Dv1_j2kLVLH = (xx_j2kLVLH_rel_transfer(1,4:6) - xx_j2kLVLH_rel_Last(end,4:6))*1e3;
Dv2_j2kLVLH = (xx_j2kLVLH_rel_Next2(1,4:6) - xx_j2kLVLH_rel_transfer(end,4:6))*1e3;
Dv12_j2kLVLH = [Dv1_j2kLVLH;Dv2_j2kLVLH];
dv12_j2kLVLH = sqrt(sum(Dv12_j2kLVLH.^2,2));

% 如果末端位置误差大于1m，则将脉冲均记为无穷
if rf_error>1e-3
    warning('优化未收敛，差%.1e km未到达目标点',norm(rf_error))
    Dv12_j2kLVLH = Inf*Dv12_j2kLVLH;
    dv12_j2kLVLH = [Inf,Inf];
end

minDist = min(sqrt(sum(xx_j2kLVLH_rel_transfer(:,1:3).^2,2)));

%% 输出
% 存自然轨道段
% xf_j2k_target_Last = T_Rot2ECJ2k(aux.jd0+tmf/86400,xf_MCR_target_Last,aux.C_Mat, 'MCEMR',0);
% xf_MCR_chaser_Last = xf_MCR_target_Last + xx_MCR_rel_Last(end,:);
% xf_j2k_chaser_Last = T_Rot2ECJ2k(aux.jd0+tmf/86400,xf_MCR_chaser_Last,aux.C_Mat, 'MCEMR',0);
% xf_j2k_rel_Last = xf_j2k_chaser_Last - xf_j2k_target_Last;
% xf_j2kLVLH_rel_Last = T_TCR2TCO_eph(xf_j2k_rel_Last,xf_j2k_target_Last,[],'LVLH');

auxSFF_dv = auxSFF;
try % 去除结构体不必要的信息
    auxSFF_dv = rmfield(auxSFF_dv, {'targetObserveDirec'});
catch
end
auxSFF_dv.naturalSFF(ii_segment).t0 = tmf;
auxSFF_dv.naturalSFF(ii_segment).x0_MCRLVLH_rel = xf_MCRLVLH_rel_Next;
auxSFF_dv.naturalSFF(ii_segment).x0_MCR_target = xx_MCR_target_Next(end,:);
auxSFF_dv.naturalSFF(ii_segment).x0_MCR_chaser = xx_MCR_chaser_Next(end,:);
auxSFF_dv.naturalSFF(ii_segment).x0_j2kLVLH_rel = xx_j2kLVLH_rel_Next2(1,:);
auxSFF_dv.naturalSFF(ii_segment).x0_j2k_target = xx_j2k_target_Next2(1,:);
auxSFF_dv.naturalSFF(ii_segment).x0_j2k_chaser = xx_j2k_chaser_Next2(1,:);
if ii_segment > 1
    auxSFF_dv.naturalSFF(ii_segment-1).tf = tm0;
end

% 存转移脉冲段
if ii_segment == 1
    auxSFF_dv.transferSFF(ii_segment).tm0 = epoch_SSF_Last;
    auxSFF_dv.transferSFF(ii_segment).tmf = tm0;
    auxSFF_dv.transferSFF(ii_segment).x0_j2k_target = xx_j2k_target_Last(1,:);
    auxSFF_dv.transferSFF(ii_segment).x0_j2k_chaser_beforeDv = xx_j2k_chaser_Last(1,:);
    dv0_j2k = T_TCO2TCR_eph(dv0,xx_j2k_target_Last(1,:),aa_j2k_target_Last(1,:),'LVLH');
    auxSFF_dv.transferSFF(ii_segment).x0_j2k_chaser_beforeDv(4:6) = xx_j2k_chaser_Last(1,4:6)-dv0_j2k*1e-3;
    auxSFF_dv.transferSFF(ii_segment).dvm0_j2kLVLH = dv0;
    auxSFF_dv.transferSFF(ii_segment).dvm0_j2k = dv0_j2k;
    auxSFF_dv.transferSFF(ii_segment).dvm0_MCRLVLH = xx_MCRLVLH_rel_Last(1,4:6);
    auxSFF_dv.transferSFF(ii_segment).dvm0_MCR = xx_MCR_rel_Last(1,4:6);
end

auxSFF_dv.transferSFF(ii_segment+1).tm0 = tm0;
auxSFF_dv.transferSFF(ii_segment+1).tmf = tmf;

Dv1_MCR = T_TCO2TCR_eph(Dv1_MCRLVLH,rv0f_MCRLVLH_rel_mane(1,:),[],'LVLH');
Dv2_MCR = T_TCO2TCR_eph(Dv2_MCRLVLH,rv0f_MCRLVLH_rel_mane(2,:),[],'LVLH');
Dv1_j2k = T_TCO2TCR_eph(Dv1_j2kLVLH,xx_j2k_target_transfer(1,:),[],'LVLH');
Dv2_j2k = T_TCO2TCR_eph(Dv1_j2kLVLH,xx_j2k_target_transfer(end,:),[],'LVLH');

auxSFF_dv.transferSFF(ii_segment+1).x0_j2k_target = xx_j2k_target_transfer(1,:);
auxSFF_dv.transferSFF(ii_segment+1).x0_j2k_chaser_beforeDv = xx_j2k_chaser_transfer(1,:);
auxSFF_dv.transferSFF(ii_segment+1).x0_j2k_chaser_beforeDv(4:6) = xx_j2k_chaser_transfer(1,4:6)-Dv1_j2k*1e-3;

auxSFF_dv.transferSFF(ii_segment+1).dvm0_j2kLVLH = Dv1_j2kLVLH;
auxSFF_dv.transferSFF(ii_segment+1).dvm0_j2k = Dv1_j2k;
auxSFF_dv.transferSFF(ii_segment+1).dvm0_MCRLVLH = Dv1_MCRLVLH;
auxSFF_dv.transferSFF(ii_segment+1).dvm0_MCR = Dv1_MCR;
auxSFF_dv.transferSFF(ii_segment+1).dvmf_j2kLVLH = Dv2_j2kLVLH;
auxSFF_dv.transferSFF(ii_segment+1).dvmf_j2k = Dv2_j2k;
auxSFF_dv.transferSFF(ii_segment+1).dvmf_MCRLVLH = Dv2_MCRLVLH;
auxSFF_dv.transferSFF(ii_segment+1).dvmf_MCR = Dv2_MCR;
auxSFF_dv.transferSFF(ii_segment+1).rf_error = rf_error;
auxSFF_dv.transferSFF(ii_segment+1).minDist = minDist;

%% 画图
if isplot == 1
    % 前后段各截取一段轨道
%     if ii_segment == 1
%         index_Last = 1:500;
%     else
%         index_Last = 1100:1500;
%     end
%     index_Next2 = 1:600;
    index_Last = 1:length(xx_MCRLVLH_rel_Last);
    index_Next2 = 1:length(xx_MCRLVLH_rel_Next2);
    %---------月心旋转系下的LVLH-----------
    figure(1)
    hold on;
    p2 = plot3(xx_MCRLVLH_rel_transfer(:,1), xx_MCRLVLH_rel_transfer(:,2), xx_MCRLVLH_rel_transfer(:,3),'Color','r');
    hold on;
    p3 = plot3(xx_MCRLVLH_rel_transfer(1,1), xx_MCRLVLH_rel_transfer(1,2), xx_MCRLVLH_rel_transfer(1,3),'g^');
    p4 = plot3(xx_MCRLVLH_rel_transfer(end,1), xx_MCRLVLH_rel_transfer(end,2), xx_MCRLVLH_rel_transfer(end,3),'rv');
    if ii_segment >= 3
    p5 = plot3(xx_MCRLVLH_rel_Next2(index_Next2,1), xx_MCRLVLH_rel_Next2(index_Next2,2), xx_MCRLVLH_rel_Next2(index_Next2,3),'Color',[237, 177, 32]/255);
    end
    p6 = plot3(0,0,0,'ks','MarkerSize',5);
    if ii_segment == 1
        p1 = plot3(xx_MCRLVLH_rel_Last(index_Last,1), xx_MCRLVLH_rel_Last(index_Last,2), xx_MCRLVLH_rel_Last(index_Last,3),'Color',[16, 63, 145]/255);
%         legend([p1,p2,p5,p3,p4],'分离轨道','转移轨道','任务轨道1','变轨起始点','变轨末端点','Location','northeastoutside');
    else
        p1 = plot3(xx_MCRLVLH_rel_Last(index_Last,1), xx_MCRLVLH_rel_Last(index_Last,2), xx_MCRLVLH_rel_Last(index_Last,3),'Color',[237, 177, 32]/255);
        legend([p1,p2,p3,p4],'任务轨道','转移轨道','变轨起始点','变轨末端点','Location','northeastoutside');
    end
    box on; grid on; grid minor; hold off; 
    axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); zlabel('\itz_L \rm[km]')
    scale_ref = naturalSFF(ii_segment).scale_ref;
    yLimit = scale_ref*((auxSFF.flag_refDirec == -1)*[-1.3,0.2] + (auxSFF.flag_refDirec == 1)*[-0.2,1.3]);
    ylim(yLimit); xlim(scale_ref*[-1,1]);
%     ylim(naturalSFF(ii_segment).scale_ref*[-1.3,0.2]); xlim(naturalSFF(ii_segment).scale_ref*[-1,1]);
    set(gca,'FontSize',15); % 设置文字大小，同时影响坐标轴标注、图例、标题等。
    title('MCR LVLH')
    view(0,90)

    %---------地心惯性系下的LVLH-----------
    figure(2)
    hold on
    p2 = plot3(xx_j2kLVLH_rel_transfer(:,1), xx_j2kLVLH_rel_transfer(:,2), xx_j2kLVLH_rel_transfer(:,3),'Color','r');
    hold on
    p3 = plot3(xx_j2kLVLH_rel_transfer(1,1), xx_j2kLVLH_rel_transfer(1,2), xx_j2kLVLH_rel_transfer(1,3),'g^');
    p4 = plot3(xx_j2kLVLH_rel_transfer(end,1), xx_j2kLVLH_rel_transfer(end,2), xx_j2kLVLH_rel_transfer(end,3),'rv');
    if ii_segment >= 3
    p5 = plot3(xx_j2kLVLH_rel_Next2(index_Next2,1), xx_j2kLVLH_rel_Next2(index_Next2,2), xx_j2kLVLH_rel_Next2(index_Next2,3),'Color',[237, 177, 32]/255);
    end
    plot3(0,0,0,'ks','MarkerSize',5); 
    if ii_segment == 1
        p1 = plot3(xx_j2kLVLH_rel_Last(index_Last,1), xx_j2kLVLH_rel_Last(index_Last,2), xx_j2kLVLH_rel_Last(index_Last,3),'Color',[16, 63, 145]/255);
%         legend([p1,p2,p5,p3,p4],'分离轨道','转移轨道','任务轨道1','变轨起始点','变轨末端点','Location','northeastoutside');
    else
        p1 = plot3(xx_j2kLVLH_rel_Last(index_Last,1), xx_j2kLVLH_rel_Last(index_Last,2), xx_j2kLVLH_rel_Last(index_Last,3),'Color',[237, 177, 32]/255);
%         legend([p1,p2,p5,p3,p4],['任务轨道',num2str(ii_segment-1)],['任务轨道',num2str(ii_segment)],'转移轨道','变轨起始点','变轨末端点','Location','northeastoutside');
        legend([p1,p2,p3,p4],'任务轨道','转移轨道','变轨起始点','变轨末端点','Location','northeastoutside');
    end
    box on; grid on; grid minor; hold off; 
    axis equal; xlabel('\itx_L \rm[km]'); ylabel('\ity_L \rm[km]'); zlabel('\itz_L \rm[km]')
    xlim(1.2*[min(xx_j2kLVLH_rel_Next2(:,1)),max(xx_j2kLVLH_rel_Next2(:,1))]); 
    ylim(1.1*[min(xx_j2kLVLH_rel_Next2(:,2)),max(xx_j2kLVLH_rel_Next2(:,2))]); 
    zlim(1.1*[min(xx_j2kLVLH_rel_Next2(:,2)),max(xx_j2kLVLH_rel_Next2(:,2))]); 
    title('ECI LVLH')
    set(gca,'FontSize',15);
    view(0,90)
%     error_rf_MCRLVLH = xx_MCRLVLH_rel_Next2(1,1:3)-xx_MCRLVLH_rel_transfer(end,1:3);
%     error_rf_j2kLVLH = xx_j2kLVLH_rel_Next2(1,1:3)-xx_j2kLVLH_rel_transfer(end,1:3);
end
