dbstop if error

clear
% close all
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg; format compact; 
% warning off
set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字
set(0,'defaultLineLineWidth',1.5)

CoreNum = 63; %设定机器CPU核心数量
if isempty(gcp('nocreate')) %如果并行未开启
    parpool(CoreNum);
end

%% 加载任务轨道数据
isplot = 0;
load('NaturalSSF_DRO3_HighFidelity_MCR1')

%% 转移轨道段1：从DRO分离，再转移至短距离相对运动任务轨道
% DRO to 任务1: time < 1 days, impulse < 10 m/s, seperation impulse 0.5-1 m/s
ii_segment = 1;
epoch_sep = auxSFF.naturalSFF(ii_segment).t0;
x0_MCR_target = auxSFF.naturalSFF(ii_segment).x0_MCR_target;

% 为保证安全，计算变轨容许的分离脉冲方向
% dv0_theta = 90; % [0,180] deg,面外角度，为减小脉冲，最好设成90°,此时z向分量为0
% dv0_phi = mean(phi_bounds); % phi_bounds deg,面内角度，与分离时刻的主星位置有关

% 循环计算距离和脉冲大小
% phi_all = linspace(phi_bounds(1)-90,phi_bounds(2)+90,100);
% dv1_j2kLVLH_all = zeros(100,2);
% minDist_all = zeros(100,1);
% isplot = 0;
% parfor ii_phi = 1:100
% dv0_phi = phi_all(ii_phi);

dv0_norm = 0.5; % [0.5,1] m/s
dv0 = dv0_norm * auxSFF.targetObserveDirec; % m/s
% dv_drt = [1,1,0];
% dv0 = dv0_norm * dv_drt/norm(dv_drt); % m/s
dt_coast = 86400*0.1; % 分离后滑行时间,sec
% 变轨设置
dt_transfer = 86400*1; % 轨道转移时间
tm10 = epoch_sep + dt_coast;
tm1f = tm10 + dt_transfer;
% dt_transfer = tmf-tm0; < 1 day
tspan_transfer1 = [tm10,tm1f]; % 分离后滑行时间
auxSFF1 = NaturalSFFTransfer(auxSFF,tspan_transfer1,ii_segment,dv0,isplot,DynamicModel);


%% 转移轨道段2：从短距离任务轨道转移至中距离任务轨道
% 任务1 to 任务2: time < 5 days, impulse < 10 m/s
ii_segment = 2;
dt_transfer2 = 1*86400;
epoch_segment2 = auxSFF1.naturalSFF(ii_segment).t0;
tm20 = epoch_segment2 - dt_transfer2/2;
tm2f = epoch_segment2 + dt_transfer2/2;
% dt_transfer = tmf-tm0;
tspan_transfer2 = [tm20,tm2f]; % 分离后滑行时间

auxSFF2 = NaturalSFFTransfer(auxSFF1,tspan_transfer2,ii_segment,[],isplot,DynamicModel);


%% 转移轨道段3：从中距离任务轨道转移至远距离任务轨道
% 任务2 to 任务3: time < 10 days, impulse < 20 m/s
ii_segment = 3;
dt_transfer3 = 3*86400; % sec
epoch_segment3 = auxSFF2.naturalSFF(ii_segment).t0;
tm30 = epoch_segment3 - dt_transfer3/2;
tm3f = epoch_segment3 + dt_transfer3/2;
% dt_transfer = tmf-tm0;
tspan_transfer3 = [tm30,tm3f]; % 分离后滑行时间

auxSFF_final = NaturalSFFTransfer(auxSFF2,tspan_transfer3,ii_segment,[],isplot,DynamicModel);

%% 存储数据
aux_SFF = rmfield(auxSFF_final, {'flag_refDirec'});
aux_SFF.naturalSFF = rmfield(aux_SFF.naturalSFF, {'scale_ref','MaxDist'});
aux_SFF.transferSFF = rmfield(aux_SFF.transferSFF, {'rf_error','minDist'});
aux_DynModel = [];
aux_DynModel.jd0 = aux_SFF.jd0;
aux_DynModel.t0UTC = [];
if strcmp(DynamicModel,'HighFidelity')
    aux_DynModel = SPICEinitialize(aux_DynModel,1); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux_DynModel = initialize(aux_DynModel); % 初始化
end
aux_SFF.jd0UTC = juliandate(aux_DynModel.t0UTC);
% aux_SFF.jd0TDB = juliandate(aux_DynModel.t0TDB);
aux_DynModel = rmfield(aux_DynModel, {'jd0','t0UTC','t0TDB','jd0_2et'});
aux_SFF = rmfield(aux_SFF, {'jd0'});

for ii_segment = 1:3
    aux_SFF.naturalSFF(ii_segment).x0_j2k_A = aux_SFF.naturalSFF(ii_segment).x0_j2k_target;
    aux_SFF.naturalSFF(ii_segment).x0_j2k_B = aux_SFF.naturalSFF(ii_segment).x0_j2k_chaser;
end
for ii_segment = 1:4
    aux_SFF.transferSFF(ii_segment).x0_j2k_A = aux_SFF.transferSFF(ii_segment).x0_j2k_target;
    aux_SFF.transferSFF(ii_segment).x0_j2k_B_beforeDv = aux_SFF.transferSFF(ii_segment).x0_j2k_chaser_beforeDv;
end
aux_SFF.x00_j2k_A = aux_SFF.x00_j2k_target;

aux_SFF = rmfield(aux_SFF, {'x00_j2k_target'});
aux_SFF.naturalSFF = rmfield(aux_SFF.naturalSFF, {'x0_j2kLVLH_rel','x0_MCR_chaser','x0_MCR_target','x0_MCRLVLH_rel','x0_j2k_chaser','x0_j2k_target'});
aux_SFF.transferSFF = rmfield(aux_SFF.transferSFF, {'dvm0_j2kLVLH','dvm0_MCR','dvm0_MCRLVLH','dvmf_j2kLVLH','dvmf_MCR','dvmf_MCRLVLH','x0_j2k_target','x0_j2k_chaser_beforeDv'});

% save data_DROSFF aux_SFF aux_DynModel

