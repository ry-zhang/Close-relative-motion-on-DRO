%% DRO近距离绕飞编队优化程序（MATLAB-STK联调）
% v0：杨驰航(2022/05/16)，email: ychhtl@foxmail.com
%     创建并撰写核心程序

%% 程序说明
% auxSFF = DROSFF_OrbitOptimize(auxSFF,DynamicModel)
% 
% 输入参数:
% -------------------------------------------------------------------------
% auxSFF.jd0:                           DRO任务入轨时刻(TBD儒略日)
%       .x00_j2k_AB:                    DRO任务入轨时刻轨道状态
%       .coasting:                      储存滑行段轨道数据的结构体
%       .naturalSFF:                    储存三次编队轨道数据的结构体
%       .naturalSFF(i).scale_ref:       第i个编队的轨道尺度
%       .naturalSFF(i).t0:              第i个编队的优化起始时刻(相对于jd0的秒数)
%       .naturalSFF(i).tf:              第i个编队的优化终端时刻(相对于jd0的秒数)
%       .naturalSFF(i).x0_MCEMR_A       第i个编队的起始时刻A星状态(月心地月旋转系)
%       .naturalSFF(i).x0_j2k_A         第i个编队的起始时刻A星状态(地心J2000坐标系)
%       .naturalSFF(i).t_sample         第i个编队的轨道采样时刻(地心J2000坐标系)
%       .naturalSFF(i).xx_MCEMR_A       第i个编队的A星轨道数据(月心地月旋转系)
%       .naturalSFF(i).xx_j2k_A         第i个编队的A星轨道数据(地心J2000坐标系)
% 
% auxSTK.scenario:                      STK的场景接口
%       .sat_AB_coasting                STK的滑行段卫星接口
%       .sat_forma                      STK的编队任务段卫星接口
% 
% 输出参数:
% -------------------------------------------------------------------------
% auxSFF.(...)                          除输入状态外加入优化结果
%       .naturalSFF(i).x0_j2k_B         第i个编队的起始时刻B星状态(月心地月旋转系)
%       .naturalSFF(i).x0_MCEMR_B       第i个编队的起始时刻B星状态(地心J2000坐标系)
%       .naturalSFF(i).xx_j2k_B         第i个编队的B星轨道数据(月心地月旋转系)
%       .naturalSFF(i).xx_MCEMR_B       第i个编队的B星轨道数据(地心J2000坐标系)
%       .naturalSFF(i).x0_j2kLVLH_REL   第i个编队的起始时刻B星状态(地心J2000LVLH)
%       .naturalSFF(i).x0_MCEMRLVLH_REL 第i个编队的起始时刻B星状态(月心地月旋转系LVLH)
%       .naturalSFF(i).xx_j2kLVLH_REL   第i个编队的相对运动轨道数据(地心J2000LVLH)
%       .naturalSFF(i).xx_MCEMRLVLH_REL 第i个编队的相对运动轨道数据(月心地月旋转系LVLH)
%       .naturalSFF(i).MaxDist          第i个编队的相对运动轨迹距CRTBP下参考轨道的最大距离
% 
% 程序中航天器状态变量命名规则
% -------------------------------------------------------------------------
% 命名格式: 变量名_坐标系_航天器_后缀
% 变量名:   x0-初始状态，xf-末端状态，xx-整个轨道状态，对于a/r/v，采用同样命名规则
% 坐标系:   MCEMR-月心地月旋转系，j2k-地心J2000坐标系，MCEMRLVLH-MCEMR下的LVLH相对系，j2kLVLH-j2k下的LVLH相对系
% 航天器:   A-A星，B-B星，REL-两者相对运动
% 后缀:     可用于标识额外信息，可有多段。例如描述不同阶段：Last\Next\transfer,于描述不同模型：CRTBP/eph
% 例:       xx_MCEMR_A_Next \ x0_MCEMRLVLH_REL_Last \ xf_j2kLVLH_REL_Last

%% 主程序
function auxSFF = DROSFF_OrbitOptimize(auxSFF,auxSTK)

naturalSFF = auxSFF.naturalSFF;

% ---------------------------------设定参考轨道选取方向-----------------------
flag_refDirec = 1;

%% 计算CRTBP下的参考轨道
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% --------------------------------加载CRTBP下的轨道数据-------------------------
load('DRO1to2_FloquetEigen')
x0_DRO_loop = x0_DRO_M_3d;
x0_REL_loop = flag_refDirec/con.r_norma/EigenVector.p3(2)*EigenVector.p3';

% -------------------------------------CRTBP下积分--------------------------------
tf = para.T0;
sol_CRTBP = ode113(@(t,x)eom_rel3b(t,x,con.mu),[0 tf], [x0_DRO_loop, x0_REL_loop], opts_ode);
t_sample_CRTBP = linspace(0,tf,2000);
sol_deval_CRTBP = deval(sol_CRTBP,t_sample_CRTBP)';
xx_MCEMR_A_CRTBP = sol_deval_CRTBP(:,1:6);
rr_MCEMRLVLH_REL_CRTBP = sol_deval_CRTBP(:,7:9);
vv_MCEMRLVLH_REL_CRTBP = sol_deval_CRTBP(:,10:12);

theta_all = mod(atan2(-xx_MCEMR_A_CRTBP(:,2),xx_MCEMR_A_CRTBP(:,1)),2*pi);
rr_MCEMRLVLH_REL_CRTBP_km = rr_MCEMRLVLH_REL_CRTBP.*con.r_norma;
vv_MCEMRLVLH_REL_CRTBP_km = vv_MCEMRLVLH_REL_CRTBP.*con.v_norma;

%% DE优化星历下的自然编队构型

% 设置DE优化参数
opts_DE = struct('Max_nfes',3000,'Display',1,'Resize',0,'UseParallel',0,'PopSize',63);

% ---------------------------依次优化多条轨道-----------------------------------
num_forma = length(naturalSFF);
for ii_forma = 1:num_forma

    % 
    scale_ref = naturalSFF(ii_forma).scale_ref;

    % ----------------------------------主星轨道递推-------------------------------
    % 当前优化时段的状态
    xx_MCEMR_A_STK = naturalSFF(ii_forma).xx_MCEMR_A;

    % ----------------------------------计算参考轨道-------------------------------
    theta_MCEMR_A_STK = mod(atan2(-xx_MCEMR_A_STK(:,2),xx_MCEMR_A_STK(:,1)),2*pi);
    x_MCEMR_interp = scale_ref*interp1(theta_all,rr_MCEMRLVLH_REL_CRTBP_km(:,1),theta_MCEMR_A_STK);
    y_MCEMR_interp = scale_ref*interp1(theta_all,rr_MCEMRLVLH_REL_CRTBP_km(:,2),theta_MCEMR_A_STK);
    vx_MCEMR_interp = scale_ref*interp1(theta_all,vv_MCEMRLVLH_REL_CRTBP_km(:,1),theta_MCEMR_A_STK);
    vy_MCEMR_interp = scale_ref*interp1(theta_all,vv_MCEMRLVLH_REL_CRTBP_km(:,2),theta_MCEMR_A_STK);
    rr_MCEMRLVLH_REL_ref = [x_MCEMR_interp,y_MCEMR_interp,zeros(size(x_MCEMR_interp))];

    % ----------------------------------设置优化上下界------------------------------
    r_min = scale_ref*min(rr_MCEMRLVLH_REL_CRTBP_km(:,[1,2])); 
    r_max = scale_ref*max(rr_MCEMRLVLH_REL_CRTBP_km(:,[1,2]));
    v_min = scale_ref*min(vv_MCEMRLVLH_REL_CRTBP_km(:,[1,2])); 
    v_max = scale_ref*max(vv_MCEMRLVLH_REL_CRTBP_km(:,[1,2]));
    scale_factor = [1,1,1e5,1e5]; % 优化归一化参数
    r0 = [x_MCEMR_interp(1),y_MCEMR_interp(1)];
    v0 = [vx_MCEMR_interp(1),vy_MCEMR_interp(1)];
    lb = scale_factor.*[r0-(r_max-r_min)/4, v0-(v_max-v_min)/4];
    ub = scale_factor.*[r0+(r_max-r_min)/4, v0+(v_max-v_min)/4];

    % -----------------------------------设置优化函数-----------------------------------
    satA = auxSTK.sat_forma(ii_forma).satA;
    satB = auxSTK.sat_forma(ii_forma).satB;
    fx = @(x)DROforma_dist(x,rr_MCEMRLVLH_REL_ref,scale_factor,satB,satA);

    % --------------------------------------DE优化--------------------------------------
    tic
    [MaxDist_history,x0_MCEMRLVLH_scaled_history] = PfhDE(fx,lb,ub,opts_DE);
    t_end = toc;
    
    % --------------------------计算并存储B星与相对运动的初值与轨道---------------------
    t0 = naturalSFF(ii_forma).t0;
    tf = naturalSFF(ii_forma).tf;
    tstep = 3600;
    naturalSFF_ii = SaveOrbit2Aux(naturalSFF(ii_forma),satA.sat,satB.sat,t0,tf,tstep);
    filedname = fieldnames(naturalSFF_ii);
    for ii_k = 1:length(filedname)
        eval(['naturalSFF(ii_forma).',filedname{ii_k},' = naturalSFF_ii.',filedname{ii_k},';'])
    end
    
    naturalSFF(ii_forma).MaxDist = MaxDist_history(end);
    
    % --------------------------------显示优化结果-------------------------------
    if opts_DE.Display == 1
        disp(['   优化次数' , '    最大偏差 / 轨道尺度', '      所耗时间']);
    else
        if ii_forma == 1
            disp(['   优化次数' , '    最大偏差 / 轨道尺度', '      所耗时间']);
        end
    end
    disp(['      ',num2str(ii_forma), '           ',num2str(MaxDist_history(end),'%5.2f'),...
        ' / ',num2str(mean(scale_ref*0.7)),' km        ', num2str(t_end,'%04.2f'),' s']);
    
end
auxSFF.naturalSFF = naturalSFF;
end