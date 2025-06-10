% 根据mainSFFtransfer输出的aux_SFF，计算整个编队飞行过程中的数据

dbstop if error

clear
close all
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')
load('data_DROSFF.mat')

%% 创建aux
aux = []; 
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.t0UTC = datetime(aux_SFF.jd0UTC,'ConvertFrom','juliandate','Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
DynamicModel = 'HighFidelity';
if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,0); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux = initialize(aux); % 初始化
else
    error('Wrong DynamicModel')
end

naturalSFF = aux_SFF.naturalSFF;
transferSFF = aux_SFF.transferSFF;

figure(1)
hold on; axis equal
%% 滑行段
tspan_sec = [0,transferSFF(1).tm0];
% t_sample = linspace(tspan_sec(1),tspan_sec(2),max(500,ceil(diff(tspan_sec)/60)+1));
t_sample = tspan_sec(1):60:tspan_sec(2);
x0_j2k_A = aux_SFF.x00_j2k_A;
[xx_j2k_A,~] = Propagate_Ephj2k(x0_j2k_A,tspan_sec,t_sample,aux);
aux_SFF.Coasting.t_sample = t_sample;
aux_SFF.Coasting.xx_j2k_AB = xx_j2k_A;

%% 自然段
for ii_NatSegment = 1:3
    tspan_sec = [naturalSFF(ii_NatSegment).t0,naturalSFF(ii_NatSegment).tf];
    t_sample = tspan_sec(1):60:tspan_sec(2);
    
    % 主星
    x0_j2k_A = naturalSFF(ii_NatSegment).x0_j2k_A;
    [xx_j2k_A,~] = Propagate_Ephj2k(x0_j2k_A,tspan_sec,t_sample,aux);
    naturalSFF(ii_NatSegment).t_sample = t_sample;
    naturalSFF(ii_NatSegment).xx_j2k_A = xx_j2k_A;
    
    % 副星
    x0_j2k_B = naturalSFF(ii_NatSegment).x0_j2k_B;
    [xx_j2k_B,~] = Propagate_Ephj2k(x0_j2k_B,tspan_sec,t_sample,aux);
    naturalSFF(ii_NatSegment).xx_j2k_B = xx_j2k_B;
    naturalSFF(ii_NatSegment).xx_j2k_rel = xx_j2k_B-xx_j2k_A;
    
    figure(1)
    hold on
    plot3(naturalSFF(ii_NatSegment).xx_j2k_rel(:,1),naturalSFF(ii_NatSegment).xx_j2k_rel(:,2),naturalSFF(ii_NatSegment).xx_j2k_rel(:,3),'b')
end

%% 变轨段
for ii_TranSegment = 1:4
    
    tspan_sec = [transferSFF(ii_TranSegment).tm0,transferSFF(ii_TranSegment).tmf];
    t_sample = tspan_sec(1):60:tspan_sec(2);
    
    % 主星
    x0_j2k_A = transferSFF(ii_TranSegment).x0_j2k_A;
    x0_j2k_B_beforeDv = transferSFF(ii_TranSegment).x0_j2k_B_beforeDv;
    [xx_j2k_A,~] = Propagate_Ephj2k(x0_j2k_A,tspan_sec,t_sample,aux);
    transferSFF(ii_TranSegment).t_sample = t_sample;
    transferSFF(ii_TranSegment).xx_j2k_A = xx_j2k_A;
    
    % 副星
    x0_j2k_B_afterDv = [x0_j2k_B_beforeDv(1:3), x0_j2k_B_beforeDv(4:6)+transferSFF(ii_TranSegment).dvm0_j2k*1e-3];
    [xx_j2k_B,~] = Propagate_Ephj2k(x0_j2k_B_afterDv,tspan_sec,t_sample,aux);
    transferSFF(ii_TranSegment).xx_j2k_B = xx_j2k_B;
    transferSFF(ii_TranSegment).xx_j2k_rel = xx_j2k_B-xx_j2k_A;
    
    figure(1)
    plot3(transferSFF(ii_TranSegment).xx_j2k_rel(:,1),transferSFF(ii_TranSegment).xx_j2k_rel(:,2),transferSFF(ii_TranSegment).xx_j2k_rel(:,3),'r')
end

aux_SFF.naturalSFF = naturalSFF;
aux_SFF.transferSFF = transferSFF;

save data_DROSFF_withOrbit2 aux_SFF aux_DynModel