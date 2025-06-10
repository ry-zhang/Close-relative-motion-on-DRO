% ------------------------------编队任务说明-----------------------------
% 整个编队任务过程分为八个轨道段，其中包括滑行段、DRO分离段、三个任务轨道段、以及三个转移轨道段
% 按时序依次为: 分离段 - 转移1 - 任务1 - 转移2 - 任务2 - 转移3 - 任务3
% 滑行段：从初值滑行至分离窗口（分离时刻）
% 分离段：分离末端点 to 任务1: time < 1 days
% 任务1:               duration > 30 days, distance [5,50] km
% 任务1 to 任务2:      time < 5 days
% 任务2:               duration > 60 days, distance [50,100]km
% 任务2 to 任务3:      time < 10 days
% 任务3:               duration > 90 days, distance > 100 km

% 本部分首先优化计算三个任务轨道段，其次优化计算中间的轨道转移段

dbstop if error

clear
close all
addpath('./subFunc(CRTBP)')
addpath('./subFunc(STK)')
addpath('./subFunc(optimize)')
format longg; format compact

set(0,'defaultAxesFontName', 'TimesSimSun','defaultTextFontName', 'TimesSimSun');
set(0,'defaultAxesFontSize',15,'defaultTextFontSize',15)
set(0,'defaultLineLineWidth',1.5)

% CoreNum = 63; %设定机器CPU核心数量
% if isempty(gcp('nocreate')) %如果并行未开启
%     parpool(CoreNum);
% end
T_Moon = (23*60+56-52.6)*60; %月球上中天的大致周期
%% 设定编队参数
% ------------------------------------设置初始历元,星历积分初值------------------------
% auxSFF.jd0 = 2460109.81760632; % TDB
% auxSFF.x00_j2k_AB = [260861.357166097,167404.871364748,78987.9106263715,...
%     -0.797706314601184,0.975263733684134,0.530154969412835]; % 地心J2K

auxSFF.jd0 = 2460369.15157369; % TDB
auxSFF.x00_j2k_AB = [-303544.120404064422473, -137517.1191382767283358, -74534.7461937523330562,...
    0.5177741054759463, -1.0350865118989347, -0.5917614219330256]; % 地心J2K

dt1 = 32*T_Moon;  % 1d裕量 + 30d任务
dt2 = 66*T_Moon;  % 5d裕量 + 60d任务
% dt3 = 100*86400; % 10d裕量 + 90d任务

% ------------------------设置每一段编队飞行的尺度，尺度随scale_ref线性放缩------------
% ----------------------scale_ref = 1对应于参考轨道星间距离约[0.67,1] km--------------
s1 = 40; % 5-50km   近距离编队
s2 = 90; % 50-100km 中距离编队
% s3 = 1000; % >500km 远距离编队

% j2kLVLH中的目标相位约束
alpha_tar = [40,220]; % deg

% ------------------------------------设置轨道转移的时间------------------------
dv_sep = 0.66; % 分离脉冲大小，m/s
dt_sep = 4/24*86400; % 分离脉冲施加后自由漂移时间
dt1_transfer = 1*T_Moon; % 分离末端至任务1转移时间
dt2_transfer = 2*T_Moon; % 任务1末端至任务2转移时间
% dt3_transfer = 3*86400; % 任务2末端至任务3转移时间

% ---------------------------------将输入参数整合为向量，方便引用----------------------
% dt_all_0 = [dt1,dt2,dt3];
% s_ref = [s1,s2,s3];
% dt_transfer = [dt1_transfer,dt2_transfer,dt3_transfer];
dt_all_0 = [dt1,dt2];
s_ref = [s1,s2];
dt_transfer = [dt1_transfer,dt2_transfer];

%% 建立STK场景
% Launch a new instance of STK11 and grab it
stkMode = 1;
if stkMode == 1
    % 方法1，可以控制stk程序是否可见,调试方便
    uiapp = actxserver('STK11.application');
    uiapp.Visible = 1;
%     uiapp.Windows.Item('3D Graphics 1 - Earth').Close;
%     uiapp.Windows.Item('2D Graphics 1 - Earth').Close;
    root = uiapp.Personality2;
elseif stkMode == 2
    % 方法2，stk在后台运行,速度快
    STKXApplication = actxserver('STKX11.application');
    root = actxserver('AgStkObjects11.AgStkObjectRoot');
end

% create a new auxSTK.scenario
auxSTK.scenario = root.Children.New('eScenario','DROFormation');
% set the time period
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('JED'); 
auxSTK.scenario.SetTimePeriod(num2str(auxSFF.jd0), ['+', num2str(365), ' days']);
% reset the animation
root.ExecuteCommand('Animate * Reset');
% set the unit of time as Epoch second to convinient the calculation
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec'); 
% 创建月心旋转坐标系
CreateMCRframe(auxSTK.scenario,root)
auxSFF.PropagatorName = 'CSU_HPOP_test';
CreateProp(auxSTK.scenario,auxSFF)

%% 编队优化
% 定义alpha为B星在ECILVLH中的相位，从ECILVLH的x轴正向逆时针旋转的角度，范围为0~360
% 定义theta为A星在MCR中的相位，从MCR的x轴正向逆时针旋转的角度，范围为0~360
% 定义为逆时针，是因为MCR中的DRO与ECILVLH中的DRO周期相对运动均逆时针旋转

% 在CRTBP模型下，可寻找到theta与alpha之间的一一对应关系alpha2theta，以作为在高精度模型下的优化参考
% 优化结果表明，在高精度模型下，相同alpha对应的theta与CRTBP下相似，误差基本位于5°以内
% 因此可用来选择每一段的变轨起始点，从而确定优化时间段

% 给定两个ECILVLH中的备选分离角度alpha_tar,计算二者与入轨时刻的大致时间距离
% 考虑A星在对地定向时，在ECILVLH中的(alpha)可观测范围是30~150，与210~360
% ---------------------计算分离脉冲施加时间--------------------
% 计算入轨之后抵达某一个alpha_tarAll的最短滑行时间，建立sat_coasting
[auxSFF,auxSTK] = sat_Coasting(auxSFF,auxSTK,alpha_tar);

% --------计算入轨之后抵达theta_tar的滑行时间，在STK中建立编队卫星satA_form与satB_form-------
[auxSFF,auxSTK] = sat_naturalSFF(auxSFF,auxSTK,dt_all_0);

% --------------------------参考轨道尺寸存储至结构体---------------------
for ii_forma = 1:length(s_ref)
    auxSFF.naturalSFF(ii_forma).scale_ref = s_ref(ii_forma);
end

% ------------------建立月心地月旋转系与地心J2K下的LVLH坐标系-------------
CreateMCRLVLHframe(auxSTK.sat_forma);
CreateJ2KLVLHframe(auxSTK.sat_forma);

% -----------------------------三段任务轨道优化--------------------------
auxSFF = DROSFF_OrbitOptimize(auxSFF,auxSTK);

% ---------------根据任务轨道优化结果与转移时间，在STK中建立转移卫星,并计算滑行段----------------
[auxSFF,auxSTK] = sat_transferSFF(auxSFF,auxSTK,dt_sep,dv_sep,dt_transfer);

% -------------------------------三段转移轨道优化--------------------------
auxSFF = DROSFF_OrbitTransfer(auxSFF,auxSTK);

% -----------------------------------保存数据------------------------------
save('DRO_SSF_coasting','auxSFF')


