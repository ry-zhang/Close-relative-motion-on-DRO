%% 高精度模型下相对运动变轨计算程序
% v0：杨驰航(2022/02/22)，email: ychhtl@foxmail.com
%     创建并撰写核心程序

%% 程序说明
% auxSFF = NaturalSFFOptimize(auxSFF,DynamicModel)
% 
% Input arguments:
% -------------------------------------------------------------------------
% auxSFF.jd0:                           DRO任务入轨时刻(TBD儒略日)
%       .x00_j2k_target:                DRO任务入轨时刻轨道状态
%       .targetObserveDirec             主星分离时刻(也即naturalSFF(1).t0)在j2kLVLH中的观测指向矢量
%       .naturalSFF:                    储存三次编队轨道数据的结构体
%       .naturalSFF(i).scale_ref:       第i个编队的轨道尺度
%       .naturalSFF(i).t0:              第i个编队的优化起始时刻(相对于jd0的秒数)
%       .naturalSFF(i).tf:              第i个编队的优化终端时刻(相对于jd0的秒数)
% DynamicModel:                         动力学模型('SEMephemeris'或'HighFidelity')
% 
% Output arguments:
% -------------------------------------------------------------------------
% auxSFF.(...)                          除输入状态外加入优化结果
%       .flag_refDirec                  参考轨道选取方向
%       .naturalSFF(i).x0_MCRLVLH_rel:  t0处的MCRLVLH下的相对运动状态
%       .naturalSFF(i).x0_MCR_target:   t0处的MCR下的主星状态
%       .naturalSFF(i).x0_j2kLVLH_rel:  t0处的j2kLVLH下的相对运动状态
%       .naturalSFF(i).x0_j2k_target:   t0处的j2k下的主星状态
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
function auxSFF = NaturalSFFOptimize(auxSFF,DynamicModel)

naturalSFF = auxSFF.naturalSFF;
x00_j2k_target = auxSFF.x00_j2k_target;
jd0 = auxSFF.jd0;

% --------------------加载星历-------------------------
aux = []; 
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.jd0 = jd0;

num_segment = length(naturalSFF);
if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,1); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux = initialize(aux); % 初始化
else
    error('Wrong DynamicModel')
end

% -----------------计算第一段轨道的优化时段起始时刻的状态-----------------------
a00_j2k_target = eomj2kMtx(x00_j2k_target,0,aux);
x00_MCR_target = T_ECJ2k2Rot(aux.jd0, x00_j2k_target, a00_j2k_target,aux.C_Mat, 'MCEMR');
epoch0 = naturalSFF(1).t0;
% 由任务起始时刻，递推至第一段轨道的优化时段起始时刻(也即epoch0)
x0_MCR_target_epoch0 = Propagate_EphRotFrame(x00_MCR_target,[0,epoch0],epoch0,aux);
% ----------------------确定参考轨道选取方向-----------------------
% phi_nagetive = rDir_j2kLVLH(epoch0,x0_MCR_target_epoch0,auxSFF);
% phi_sep = atan2d(auxSFF.targetObserveDirec(2),auxSFF.targetObserveDirec(1));
% phi_diff = min(abs(phi_sep-phi_nagetive+[-360,0,360]));
% auxSFF.flag_refDirec = -1*(phi_diff<90)+1*(phi_diff>=90);

auxSFF.flag_refDirec = 1;

%% MCR下的DRO前后跟飞编队，也即CRTBP下的MCR DRO周期相对运动
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% --------------------加载CRTBP下的轨道数据-------------------------
load('FloquetEig12')
% load('FloquetEig_all.mat'); ii_DRO = 71;
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = auxSFF.flag_refDirec/con.r_norma/EigenVector.p3(2)*EigenVector.p3';%周期解

% -------------------------CRTBP下积分--------------------------------
tf = J_period_all(4,2)*2*pi;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
xx_MCR_target_CRTBP = sol_deval_CRTBP(:,1:6);
% xx_MCR_target_CRTBP(:,[1,2,4,5]) = -sol_deval_CRTBP(:,[1,2,4,5]);
rr_MCRLVLH_rel_CRTBP = sol_deval_CRTBP(:,7:9);
vv_MCRLVLH_rel_CRTBP = sol_deval_CRTBP(:,10:12);

%% 拟合CRTBP下的参考周期轨道
% ----------------------------拟合自变量与因变量--------------------------
theta_all = mod(atan2(-xx_MCR_target_CRTBP(:,2),xx_MCR_target_CRTBP(:,1)),2*pi);%计算主星在 MCR 坐标系中的几何相位。
rr_MCRLVLH_rel_CRTBP_km = rr_MCRLVLH_rel_CRTBP.*con.r_norma;
vv_MCRLVLH_rel_CRTBP_km = vv_MCRLVLH_rel_CRTBP.*con.v_norma;

% -----------------------------------拟合（最小二乘）---------------------------------
ft = fittype( 'fourier8' );%使用 8 项的 Fourier 系列进行拟合。
opts_fit = fitoptions( 'Method', 'NonlinearLeastSquares','Display','Off' );
[fitresult_x, gof_x] = fit( theta_all, rr_MCRLVLH_rel_CRTBP_km(:,1), ft, opts_fit );
[fitresult_y, gof_y] = fit( theta_all, rr_MCRLVLH_rel_CRTBP_km(:,2), ft, opts_fit );
[fitresult_vx, gof_vx] = fit( theta_all, vv_MCRLVLH_rel_CRTBP_km(:,1), ft, opts_fit );
[fitresult_vy, gof_vy] = fit( theta_all, vv_MCRLVLH_rel_CRTBP_km(:,2), ft, opts_fit );

%% DE优化星历下的自然编队构型

opts_DE = struct('Max_nfes',6000,'Display',1,'Resize',0,'UseParallel',1,'PopSize',63);
r_MCRLVLH_rel_ref=zeros(500,3);
% ---------------依次优化多条轨道-----------------------
disp(['   优化次数' , '    最大偏差 / 轨道尺度', '      所耗时间']);
for ii_segment = 1:num_segment

    scale_ref = naturalSFF(ii_segment).scale_ref;%scale_ref 轨道尺度（km）
    % 优化时间段
    epoch_SSF = naturalSFF(ii_segment).t0;
    tspan_sec = [naturalSFF(ii_segment).t0,naturalSFF(ii_segment).tf];
%     t_sample = linspace(tspan_sec(1),tspan_sec(2),max(500,ceil(diff(tspan_sec)/3600/4)+1));%它表示时间区间 [tspan_sec(1), tspan_sec(2)] 内的多个时间采样点。采样点的数量至少为 500，或者根据时间区间的长度和设定的每小时采样点数自动调整。
    t_sample = linspace(tspan_sec(1),tspan_sec(2),max(21600,ceil(diff(tspan_sec)/60)+1));%它表示时间区间 [tspan_sec(1), tspan_sec(2)] 内的多个时间采样点。采样点的数量至少为 500，或者根据时间区间的长度和设定的每小时采样点数自动调整。

    % ----------------------主星轨道递推-------------------
    % 由第一段轨道的优化时段起始时刻，递推至当前优化时段起始时刻(也即epoch_SSF)
    [x0_MCR_target,a0_MCR_target] = Propagate_EphRotFrame(x0_MCR_target_epoch0,[epoch0,epoch_SSF],epoch_SSF,aux);
    % 由当前优化时段起始时刻，递推至当前优化时段末端时刻
    [xx_MCR_target,aa_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec,t_sample,aux);

    % ----------------------计算参考轨道-------------------
    theta_MCR_targe_eph = mod(atan2(-xx_MCR_target(:,2),xx_MCR_target(:,1)),2*pi);%主星在 MCR 坐标系中的方向角/相位。
    fit_x = scale_ref*feval(fitresult_x,theta_MCR_targe_eph);%基于主星的轨道递推结果，计算参考相对轨道
    fit_y = scale_ref*feval(fitresult_y,theta_MCR_targe_eph);
    r_MCRLVLH_rel_ref = [fit_x,fit_y,zeros(size(fit_x))];%计算CRTBP下相对参考轨道（转换成eph相位）。

    % ----------------------设置优化上下界------------------
    r_min = scale_ref*min(rr_MCRLVLH_rel_CRTBP_km(:,[1,2])); 
    r_max = scale_ref*max(rr_MCRLVLH_rel_CRTBP_km(:,[1,2]));
    v_min = scale_ref*min(vv_MCRLVLH_rel_CRTBP_km(:,[1,2])); 
    v_max = scale_ref*max(vv_MCRLVLH_rel_CRTBP_km(:,[1,2]));
    scale_factor = [1,1,1e5,1e5]; % 优化归一化参数
    r0 = [fit_x(1),fit_y(1)];
    v0 = scale_ref*[feval(fitresult_vx,theta_MCR_targe_eph(1)),...
        feval(fitresult_vy,theta_MCR_targe_eph(1))];
    lb = scale_factor.*[r0-(r_max-r_min)/4, v0-(v_max-v_min)/4];%【-5000km~5000km】 %【-50m/s~50m/s】
    ub = scale_factor.*[r0+(r_max-r_min)/4, v0+(v_max-v_min)/4];

    % -----------------------设置优化函数-----------------------
    fx = @(x)dist_DROformation(x,xx_MCR_target,aa_MCR_target,...
        tspan_sec,t_sample,r_MCRLVLH_rel_ref,scale_factor,aux);

    % --------------------------DE优化--------------------------
    tic
    [MaxDist_history,x0_MCRLVLH_scaled_opt_history] = PfhDE(fx,lb,ub,opts_DE);
    t_end = toc;
    
    % ---------------------保存MCR中的优化结果----------------------------
    x0_MCRLVLH_rel_xy = x0_MCRLVLH_scaled_opt_history(end,:)./scale_factor;
    x0_MCRLVLH_rel = [x0_MCRLVLH_rel_xy(1),x0_MCRLVLH_rel_xy(2),0,x0_MCRLVLH_rel_xy(3),x0_MCRLVLH_rel_xy(4),0];
    naturalSFF(ii_segment).x0_MCRLVLH_rel = x0_MCRLVLH_rel;
    naturalSFF(ii_segment).x0_MCR_target = x0_MCR_target;
    
    % -----------------将优化结果转换到j2k中并保存-----------------------
    x0_j2k_target = T_Rot2ECJ2k(aux.jd0+tspan_sec(1)/86400,x0_MCR_target,aux.C_Mat, 'MCEMR',0);
    x0_MCR_rel = T_TCO2TCR_eph(x0_MCRLVLH_rel,x0_MCR_target,a0_MCR_target,'LVLH');
    x0_MCR_chaser = x0_MCR_target + x0_MCR_rel;
    x0_j2k_chaser = T_Rot2ECJ2k(aux.jd0+tspan_sec(1)/86400,x0_MCR_chaser,aux.C_Mat, 'MCEMR',0);
    x0_j2k_rel = x0_j2k_chaser - x0_j2k_target;
    a0_j2k_target = eomj2kMtx(x0_j2k_target,tspan_sec(1),aux);
    naturalSFF(ii_segment).x0_j2kLVLH_rel = T_TCR2TCO_eph(x0_j2k_rel,x0_j2k_target,a0_j2k_target,'LVLH');
    naturalSFF(ii_segment).x0_j2k_target = x0_j2k_target;
    
    naturalSFF(ii_segment).MaxDist = MaxDist_history(end);
    
    % ------------------------显示优化结果-----------------------
    disp(['      ',num2str(ii_segment), '           ',num2str(MaxDist_history(end),'%5.2f'),...
        ' / ',num2str(mean(scale_ref*0.7)),' km        ', num2str(t_end,'%04.2f'),' s']);
    if(ii_segment>=95)
        save('r_MCRLVLH_rel_ref.mat','r_MCRLVLH_rel_ref');
    end

    
end
save('r_MCRLVLH_rel_ref.mat','r_MCRLVLH_rel_ref')
figure(1)
% plot3(r_MCRLVLH_rel_ref(1,1),r_MCRLVLH_rel_ref(1,2),r_MCRLVLH_rel_ref(1,3),'g^');hold on;
% plot3(r_MCRLVLH_rel_ref(end,1),r_MCRLVLH_rel_ref(end,2),r_MCRLVLH_rel_ref(end,3),'r*');hold on;
plot3(r_MCRLVLH_rel_ref(:,1),r_MCRLVLH_rel_ref(:,2),r_MCRLVLH_rel_ref(:,3),'b');hold on;



auxSFF.naturalSFF = naturalSFF;
end