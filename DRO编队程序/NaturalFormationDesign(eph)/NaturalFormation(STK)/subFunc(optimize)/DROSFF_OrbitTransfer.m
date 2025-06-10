%% 高精度模型下相对运动变轨计算程序
% v0：杨驰航(2022/05/17)，email: ychhtl@foxmail.com
%     创建并撰写核心程序

%% 程序说明
% auxSFF = DROSFF_OrbitTransfer(auxSFF,auxSTK)
% 
% 输入参数:
% -------------------------------------------------------------------------
% auxSFF.jd0:                           DRO任务入轨时刻(TBD儒略日)
%       .x00_j2k_AB:                    DRO任务入轨时刻轨道状态
%       .coasting:                      储存滑行段轨道数据的结构体
%       .sepSFF:                        储存分离段轨道数据的结构体
%       .sepSFF.xx_j2kLVLH_REL:         分离段相对运动轨道数据(地心J2000LVLH)
%       .naturalSFF:                    储存三次任务段轨道数据的结构体
%       .naturalSFF(i).x0_j2kLVLH_REL   第i个编队的起始时刻B星状态(地心J2000LVLH)
%       .naturalSFF(i).xx_j2kLVLH_REL   第i个编队的相对运动轨道数据(地心J2000LVLH)
%       .transferSFF:                   储存三次变轨段轨道数据的结构体
%       .transferSFF(i).t0:             第i次转移的起始时刻(相对于jd0的秒数)
%       .transferSFF(i).tf:             第i次转移的终端时刻(相对于jd0的秒数)
% 
% auxSTK.scenario:                      STK的场景接口
%       .sat_forma                      STK的编队任务段卫星接口
%       .satB_trans                     STK的变轨段卫星接口
% 
% 输出参数:
% -------------------------------------------------------------------------
% auxSFF.(...)                          
%       .transferSFF:                   添加transferSFF中的轨道数据
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
function auxSFF = DROSFF_OrbitTransfer(auxSFF,auxSTK)

opts_fsolve = optimoptions('fsolve','Display','off',...
    'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,...
    'Algorithm','levenberg-marquardt','UseParallel',false,...
    'SpecifyObjectiveGradient', false);
for ii_forma = 1:length(auxSFF.naturalSFF)
    %% 变轨优化
    if ii_forma == 1
        x0_j2kLVLH_REL = auxSFF.sepSFF.xx_j2kLVLH_REL(end,:);
    else
        x0_j2kLVLH_REL = auxSFF.naturalSFF(ii_forma-1).xx_j2kLVLH_REL(end,:);
    end
    xf_j2kLVLH_REL_tar = auxSFF.naturalSFF(ii_forma).x0_j2kLVLH_REL;

    satB = auxSTK.satB_trans(ii_forma);
    satA = auxSTK.sat_forma(ii_forma).satA;
    
    dv0_ini = [0,0,0];
    scale_factor = 1e3;
    [dv0_scale,fval] = fsolve(@(x)DROforma_trans(x,x0_j2kLVLH_REL,xf_j2kLVLH_REL_tar,scale_factor,satB,satA),dv0_ini,opts_fsolve);
    
    dv0 = dv0_scale./scale_factor; % km/s
    [~,dvf] = DROforma_trans(dv0_scale,x0_j2kLVLH_REL,xf_j2kLVLH_REL_tar,scale_factor,satB,satA);
    
    %% 变轨数据输出与保存    
    t0_transfer = auxSFF.transferSFF(ii_forma).t0;
    tf_transfer = auxSFF.transferSFF(ii_forma).tf;
    tstep = 3600;
    auxSFF.transferSFF(ii_forma).dv0_j2kLVLH = dv0;
    auxSFF.transferSFF(ii_forma).dvf_j2kLVLH = dvf;
    auxSFF.transferSFF(ii_forma).rf_error = fval;
    
    transferSFF_ii = SaveOrbit2Aux(auxSFF.transferSFF(ii_forma),satA.sat,satB.sat,t0_transfer,tf_transfer,tstep);
    filedname = fieldnames(transferSFF_ii);
    for ii_k = 1:length(filedname)
        eval(['auxSFF.transferSFF(ii_forma).',filedname{ii_k},' = transferSFF_ii.',filedname{ii_k},';'])
    end
    
end

end