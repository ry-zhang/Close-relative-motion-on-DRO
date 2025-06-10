dbstop if error

clear
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg; format compact

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

 CoreNum = 6; %设定机器CPU核心数量
if isempty(gcp('nocreate')) %如果并行未开启
    parpool(CoreNum);
end

% 既然不能找到准确的初始ECILVLH几何相位（因为优化结果不太能确定，变动约为5°左右)，
% 所以可以考虑增加优化时间跨度，在变轨的时候再通过odeevent选取需要的相位！

%% 优化值设定
% DynamicModel = 'HighFidelity';
DynamicModel = 'SEMephemeris';%zry
iDRO = 3;
FileName = ['NaturalSSF_DRO',num2str(iDRO),'_',DynamicModel,'_MCR1'];

% ---------------------选择DRO轨道-设置初始历元,星历积分初值--------------------
switch iDRO
case 0
    % DRO1
    load('DROephMultiRev2024.mat')
    auxSFF.jd0 = jd0_temp;
    auxSFF.x00_j2k_target = x0j2k';
    % 主星在t0时刻、j2kLVLH中的观测指向矢量,期望分离矢量朝这个方向
    alpha_sepDir = 120;
    auxSFF.targetObserveDirec = [cosd(alpha_sepDir),sind(alpha_sepDir),0]; 
case 1
    % DRO1 2023.6-2027.1 <<1 hr
%     t0UTC = datetime('2023-06-14 07:36:12');
    auxSFF.jd0 = 2460109.81760632; % TBD
    auxSFF.x00_j2k_target = [260861.357166097,167404.871364748,78987.9106263715,...
        -0.797706314601184,0.975263733684134,0.530154969412835];
    % 主星在t0时刻、j2kLVLH中的观测指向矢量,期望分离矢量朝这个方向
    alpha_sepDir = 120;
    auxSFF.targetObserveDirec = [cosd(alpha_sepDir),sind(alpha_sepDir),0]; 
case 2
    % DRO2 2023.6-2027.6  ~2.5 hr
%     t0UTC = datetime('2023-06-13 23:00:16');
    auxSFF.jd0 = 2460109.45931928;
    auxSFF.x00_j2k_target = [276577.520222802,144514.098050027,67426.0628386977,...
        -0.708023138517267,1.02870635679717,0.558158622414164];
    % 主星在t0时刻、j2kLVLH中的观测指向矢量,期望分离矢量朝这个方向
    alpha_sepDir = 120;
    auxSFF.targetObserveDirec = [cosd(alpha_sepDir),sind(alpha_sepDir),0]; 
case 3
    jd0_temp = 2460109.81760632; % TBD
    x0_j2k_temp = [260861.357166097,167404.871364748,78987.9106263715,...
        -0.797706314601184,0.975263733684134,0.530154969412835];
    % 入轨时刻
    auxSFF.jd0 = juliandate(datetime('2023-01-01') + minutes(1) + seconds(9.186)); % TBD
    auxSFF.x00_j2k_target = Propagate_Ephj2k_jd0f(x0_j2k_temp,jd0_temp,auxSFF.jd0,DynamicModel);
    % 递推   
    alpha_0_deg = fit_xeci2alpha(auxSFF.jd0, auxSFF.x00_j2k_target);%eci下入轨时刻对应相位                       
    % 给定两个ECILVLH备选分离角度,计算二者与入轨时刻的大致时间距离
    alpha_tar_deg = [80,260];%eci 绝对相位
    t_0 = fit_alpha2MCRtime(alpha_0_deg); %zry初始时刻对应星历时间（eci）
    t_tar = fit_alpha2MCRtime(alpha_tar_deg); %zry分离角对应相位角（mcr）
    load('FloquetEig12')
    t_diff = mod(t_tar-t_0,pi*con.T_norma_day);%计算目标分离角和初始分离角之间的时间差，并将结果限制在一个周期内，周期映射
    [t_sort,t_sortindex] = sort(t_diff);
    if t_sort(1)>1 % 如果距离入轨时刻大于1天
        alpha_sepDir = alpha_tar_deg(t_sortindex(1));
        t0_day = t_sort(1);
    else
        alpha_sepDir = alpha_tar_deg(t_sortindex(2));
        t0_day = t_sort(2);
    end
%     t0_day = 0;
    auxSFF.targetObserveDirec = [cosd(alpha_sepDir),-sind(alpha_sepDir),0];
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
% % 短、中、远距离编队的优化时间段分别为：[t1,t2],[t2,t3],[t3,t4]
t1 = 0*86400; 
t2 = t1 + 15*86400;  % 1d转移 + (30d任务 + 10d裕量)
% t2 = t1 + 41*86400;  % 1d转移 + (30d任务 + 10d裕量)
% t3 = t2 + 75*86400;  % 5d转移 + (60d任务 + 10d裕量)
% t4 = t3 + 110*86400; % 10d转移 + (90d任务 + 10d裕量)

% t1 = t0_day*86400; 
% 
% t2 = t1 + 6.825*86400;  % 半个周期后，进入下一个窗口
% t3 = t2 + 6.825*86400;  % 半个周期后，进入下一个窗口
% t4 = t3 + 6.825*86400; % 半个周期后，进入下一个窗口
naturalSFF(1).t0 = t1; 
naturalSFF(1).tf = t2;
% naturalSFF(2).t0 = t2;
% naturalSFF(2).tf = t3;  
% naturalSFF(3).t0 = t3;
% naturalSFF(3).tf = t4; 

% -------------参考轨道尺寸，scale_ref = 1对应于参考轨道星间距离约[0.67,1] km--------
naturalSFF(1).scale_ref = 20; % 5-50km   近距离编队
% naturalSFF(2).scale_ref = 90; % 50-100km 中距离编队
% naturalSFF(3).scale_ref = 1000; % >500km 远距离编队

auxSFF.naturalSFF = naturalSFF;

% 三段任务轨道优化
auxSFF = NaturalSFFOptimize(auxSFF,DynamicModel);

% 保存数据
save(FileName,'auxSFF','DynamicModel')
save('auxSFF.mat','auxSFF')


