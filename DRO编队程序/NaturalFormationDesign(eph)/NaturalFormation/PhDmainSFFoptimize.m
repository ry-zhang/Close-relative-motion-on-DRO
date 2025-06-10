dbstop if error

clear
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg; format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

if isempty(gcp('nocreate')) %如果并行未开启
    CoreNum = 63; %设定机器CPU核心数量
    parpool(CoreNum);
end

%% 优化值设定
% DynamicModel = 'HighFidelity';
DynamicModel = 'SEMephemeris';
iDRO = 0;
FileName = ['PhDNaturalSSF_DR0',num2str(iDRO),'_',DynamicModel,'1'];

% ---------------------选择DRO轨道-设置初始历元,星历积分初值--------------------
switch iDRO
case 0
    % DRO1
    load('DROephMultiRev2024.mat')
    auxSFF.jd0 = jd0;
    auxSFF.x00_j2k_target = x0j2k';
    auxSFF.targetObserveDirec = [-1,0,0]; % 主星在t0时刻、j2kLVLH中的观测指向矢量,期望分离矢量朝这个方向
case 1
    % DRO1 2023.6-2027.1 <<1 hr
%     t0UTC = datetime('2023-06-14 07:36:12');
    auxSFF.jd0 = 2460109.81760632; % TBD
    auxSFF.x00_j2k_target = [260861.357166097,167404.871364748,78987.9106263715,...
        -0.797706314601184,0.975263733684134,0.530154969412835];
    auxSFF.targetObserveDirec = [1,0,0]; % 主星在t0时刻、j2kLVLH中的观测指向矢量,期望分离矢量朝这个方向
case 2
    % DRO2 2023.6-2027.6  ~2.5 hr
%     t0UTC = datetime('2023-06-13 23:00:16');
    auxSFF.jd0 = 2460109.45931928;
    auxSFF.x00_j2k_target = [276577.520222802,144514.098050027,67426.0628386977,...
        -0.708023138517267,1.02870635679717,0.558158622414164];
    auxSFF.targetObserveDirec = [1,0,0]; % 主星在t0时刻、j2kLVLH中的观测位置矢量,期望分离矢量朝这个方向
case 3
    % DRO2 2023.6-2027.6  ~2.5 hr
    aux.t0UTC = datetime('2023-01-01 00:00:00');
    aux = SPICEinitialize(aux);
    auxSFF.jd0 = aux.jd0;
    auxSFF.x00_j2k_target = [379402.4078133, 137717.8832331, 45812.4501215,...
        -0.59030751546, 0.68353413940, 0.34559098575];
    auxSFF.targetObserveDirec = [1,0,0]; % 主星在t0时刻、j2kLVLH中的观测位置矢量,期望分离矢量朝这个方向
otherwise
    error('no exist DRO')
end


% ------------------------------编队任务说明-----------------------------
% 整个编队任务过程分为七个轨道段，其中包括从DRO分离段、三个任务轨道段、以及三个转移轨道段
% 按时序依次为: 分离段 - 转移 - 任务1 - 转移 - 任务2 - 转移 - 任务3
% 分离段:
% 分离末端点 to 任务1: time < 1 days
% 任务1:               duration > 30 days, distance [5,50] km
% 任务1 to 任务2:      time < 5 days
% 任务2:               duration > 60 days, distance [50,100]km
% 任务2 to 任务3:      time < 10 days
% 任务3:               duration > 90 days, distance > 100 km
% 
% 三段任务轨道时间最长,为设计转移轨道,在优化自然编队任务轨道时,把优化时间调长,
% 这样编队时间段可以在变轨时间调整内自由滑动

% -------------------编队优化时段(相对于DRO入轨时刻(也即jd0)的秒数)-----------
% 短、中、远距离编队的优化时间段分别为：[t1,t2],[t2,t3],[t3,t4]
t1 = 0*86400; 
t2 = t1 + 30*86400;
t3 = t2 + 100*86400; 
t4 = t3 + 400*86400; 
naturalSFF(1).t0 = t1; 
naturalSFF(1).tf = t2;
naturalSFF(2).t0 = t2;
naturalSFF(2).tf = t3;  
naturalSFF(3).t0 = t3;
naturalSFF(3).tf = t4; 

% -------------参考轨道尺寸，scale_ref = 1对应于参考轨道星间距离约[0.67,1] km--------
naturalSFF(1).scale_ref = 6; % 4-6 km   1km 级
naturalSFF(2).scale_ref = 60; % 40-60km  10km 级
naturalSFF(3).scale_ref = 600; % 400-600km  100km 级

auxSFF.naturalSFF = naturalSFF;

% 三段任务轨道优化
auxSFF = NaturalSFFOptimize(auxSFF,DynamicModel);

% 保存数据
save(FileName,'auxSFF','DynamicModel')


