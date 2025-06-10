dbstop if error

clear
% close all
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg; format compact
warning off

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

load('FloquetEig12')
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);
opts_DE = struct('Max_nfes',6000,'Display',1,'Resize',0,'UseParallel',1,'PopSize',63);

CoreNum = 63; %设定机器CPU核心数量
if isempty(gcp('nocreate')) %如果并行未开启
    parpool(CoreNum);
end

%% 优化定义值
% 优化区间
% 为设计转移轨道，可以把编队优化时间调长，这样编队时间段可以来回滑动
numDay = 100;
tspan_sec = [0,86400*numDay];
t_sample = linspace(tspan_sec(1),tspan_sec(2),1500);
epoch_SSF_all = 0; % 编队相对于DRO的初始时刻

% 参考轨道尺寸，scale_ref = 1对应于参考轨道星间距离[0.67,1] km
% scale_ref = 20; % 5-50km   近距离编队
scale_ref = 90; % 50-100km 中距离编队
% scale_ref = 1000; % >500km 远距离编队
iDRO = 1; 

% DynamicModel = 'HighFidelity';
DynamicModel = 'SEMephemeris';

FileName = ['NaturalSSF_DRO',num2str(iDRO),'_',num2str(numDay),'day_scale',num2str(scale_ref),'_',DynamicModel];

%% CR3BP下的DRO周期相对运动
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = -1/con.r_norma/Sol_linear.vec3(2)*Sol_linear.vec3';

% 周期轨道
tf = para.T0;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
rvDRO_MCR_ref = sol_deval_CRTBP(:,1:6);

rRelMot_ref_L = sol_deval_CRTBP(:,7:9);
rRelMot_ref = rRelMot_ref_L(:,[2,3,1]);
rRelMot_ref(:,[2,3]) = -rRelMot_ref(:,[2,3]);

vRelMot_ref_L = sol_deval_CRTBP(:,10:12);
vRelMot_ref = vRelMot_ref_L(:,[2,3,1]);
vRelMot_ref(:,[2,3]) = -vRelMot_ref(:,[2,3]);

rRelMot_km_L = rRelMot_ref_L.*con.r_norma;

%% 拟合参考周期轨道
rRelMotion_km = rRelMot_ref.*con.r_norma;
vRelMotion_km = vRelMot_ref.*con.v_norma;

theta_all = mod(atan2(-rvDRO_MCR_ref(:,2),rvDRO_MCR_ref(:,1)),2*pi);
% theta_all(end) = theta_all(end);
% Set up fittype and options.
ft = fittype( 'fourier8' );
opts_fit = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts_fit.Display = 'Off';
% opts_fit.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.00074919498648];

% Fit model to data.
[fitresult_x, gof_x] = fit( theta_all, rRelMotion_km(:,1), ft, opts_fit );
[fitresult_z, gof_z] = fit( theta_all, rRelMotion_km(:,3), ft, opts_fit );

[fitresult_vx, gof_vx] = fit( theta_all, vRelMotion_km(:,1), ft, opts_fit );
[fitresult_vz, gof_vz] = fit( theta_all, vRelMotion_km(:,3), ft, opts_fit );

%% DE优化星历下的自然编队构型
% 加载星历
aux = []; 
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;

% 循环前参数设置
num_opt = length(epoch_SSF_all);
x0_TCO_optimize_all = zeros(num_opt,4);
MaxDist_all = zeros(num_opt,1);

% 选择DRO轨道-设置初始历元,星历积分初值
switch iDRO
case 1
    % DRO1
    load('DROephMultiRev2024.mat')
    aux.jd0 = jd0;
    x0_j2k_target0 = x0j2k';    
case 2
    % DRO2
    aux.t0UTC  = datetime(2023, 1, 1, 0, 0, 0); % 初始历元
    x0_j2k_target0 = [380224.0673756608739495, 140818.4956395560002420, 42079.2399940182804130, ...
        -0.5874909506324951, 0.6787788399867067, 0.3426576490741002];% 星历积分初值（j2000）
case 3
    % DRO3 2023.6-2027.1 <<1 hr
    aux.t0UTC = datetime('2023-06-14 07:36:12');
    x0_j2k_target0 = [260861.357166097,167404.871364748,78987.9106263715,...
        -0.797706314601184,0.975263733684134,0.530154969412835];
case 4
    % DRO4 2023.6-2027.6  ~2.5 hr
    aux.t0UTC = datetime('2023-06-13 23:00:16');
    x0_j2k_target0 = [276577.520222802,144514.098050027,67426.0628386977,...
        -0.708023138517267,1.02870635679717,0.558158622414164];
otherwise
    error('no exist DRO')
end

if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,1); % 初始化
else
    aux = initialize(aux); % 初始化
end
t0UTC = aux.t0UTC;

% 多次优化循环
disp(['   优化次数' , '    最大误差', '     所耗时间']);
for ii_loop = 1:num_opt
    
    % 将DRO初值向后递推epoch_SSF秒
    epoch_SSF = epoch_SSF_all(ii_loop);
    tspan_sec_re0 = tspan_sec+epoch_SSF;
    t_sample_re0 = t_sample+epoch_SSF;
    
    if epoch_SSF == 0
        xx_j2k_target0 = x0_j2k_target0;
    else
        if strcmp(DynamicModel,'HighFidelity')
            [~,xx_j2k_target0] = ode113(@eom_SPICE_stm , [0,epoch_SSF], x0_j2k_target0 , opts_ode, aux);
        else
            [~,xx_j2k_target0] = ode113(@eqom_geoMEMEJ2k, [0,epoch_SSF], x0_j2k_target0 , opts_ode, aux);
        end
    end
    
    x0_j2k_target = xx_j2k_target0(end,:);
    
%     % 初始化新的aux
%     aux = aux0; 
%     aux.t0UTC = [];
%     aux.jd0 = aux0.jd0 + epoch_SSF/86400;
%     aux = initialize(aux); % 初始化
    
    % 主星星历积分
    x0_MCR_target = T_ECJ2k2Rot(aux.jd0 + epoch_SSF/86400, x0_j2k_target, [0,0,0],aux.C_Mat, 'MCEMR');
    [xx_MCR_target,a_MCR_target] = Propagate_EphRotFrame(x0_MCR_target,tspan_sec_re0,t_sample_re0,aux);
%     figure(1); plot(xx_MCR_target(:,1), xx_MCR_target(:,2),'Color',[237, 177, 32]/255,'LineWidth',1.5); axis equal

    % 计算参考轨道
    theta_DRO_eph = mod(atan2(-xx_MCR_target(:,2),xx_MCR_target(:,1)),2*pi);
    fit_x = scale_ref*feval(fitresult_x,theta_DRO_eph);
    fit_z = scale_ref*feval(fitresult_z,theta_DRO_eph);
    rRel_ref = [fit_x,fit_z];
    
    % 设置优化上下界
    r_min = scale_ref*min(rRelMotion_km(:,[1,3])); r_max = scale_ref*max(rRelMotion_km(:,[1,3]));
    v_min = scale_ref*min(vRelMotion_km(:,[1,3])); v_max = scale_ref*max(vRelMotion_km(:,[1,3]));
    scale_factor = [1,1,1e6,1e6]; % 优化归一化参数
    r0 = [fit_x(1),fit_z(1)];
    v0 = scale_ref*[feval(fitresult_vx,theta_DRO_eph(1)),...
        feval(fitresult_vz,theta_DRO_eph(1))];
    lb = scale_factor.*[r0-(r_max-r_min)/4, v0-(v_max-v_min)/4];
    ub = scale_factor.*[r0+(r_max-r_min)/4, v0+(v_max-v_min)/4];
    
    % 设置优化函数
    fx = @(x)dist_DROformation(x,xx_MCR_target,a_MCR_target,...
        tspan_sec_re0,t_sample_re0,rRel_ref,scale_factor,aux);
    
    % 优化
    tic
    [MaxDist,x0_TCO_optim_scale] = PfhDE(fx,lb,ub,opts_DE);
%     MaxDist = 0; x0_TCO_optim_scale = zeros(1,4);
    t_end = toc;
    
    % 保存数据
    x0_TCO_optimize_all(ii_loop,:) = x0_TCO_optim_scale(end,:)./scale_factor;
    MaxDist_all(ii_loop) = MaxDist(end);
    
    % 显示优化结果
    disp(['      ',num2str(ii_loop), '      ',num2str(MaxDist(end),'%5.2f'),...
        'km     ', num2str(t_end,'%04.2f'),'s']);
end

save(FileName,'MaxDist_all', 'x0_TCO_optimize_all', 'tspan_sec',...
    't_sample', 't0UTC', 'x0_j2k_target0', 'scale_ref','epoch_SSF_all','DynamicModel')

