function [auxSFF,auxSTK] = sat_transferSFF(auxSFF,auxSTK,dt0_sep,dv0_sep_norm,dt_transfer)
% 根据任务轨道优化结果与转移时间，在STK中建立转移卫星,并更新auxSFF中的任务段的轨道数据

for ii_forma = 1:length(dt_transfer)
    % 设置forma任务段的初始数据
    if ii_forma == 1
        t0_sep = auxSFF.coasting.tf;
        tf_sep = auxSFF.coasting.tf + dt0_sep;
        
        dv0_sep = [cosd(auxSFF.alpha_tar(1)),-sind(auxSFF.alpha_tar(1)),0]*dv0_sep_norm*1e-3;
        x0_j2kLVLH_sep = [zeros(1,3),dv0_sep];
        
        % 若satB_sep已存在，则删除原卫星
        if  auxSTK.scenario.Children.Contains('eSatellite','satB_sep') 
            auxSTK.scenario.Children.Unload('eSatellite','satB_sep');
        end
        % 新建satB_sep卫星
        CoordSystem = ['Satellite/sat_A_form',num2str(ii_forma),' J2KLVLH_SC'];
        auxSTK.satB_sep = sat_initialize(auxSTK.scenario,x0_j2kLVLH_sep,CoordSystem,t0_sep,tf_sep,'satB_sep',auxSFF.PropagatorName);
        
        % 获取transfer卫星的时间跨度与初值
        t0_transfer = tf_sep;
        tf_transfer = t0_transfer + dt_transfer(1);
        x0_j2k_transfer = PosVel_j2k(auxSTK.satB_sep.sat,t0_transfer,t0_transfer,3600);
        x0_j2k_transfer = x0_j2k_transfer(2:end);
    else
        % 获取transfer卫星的时间跨度与初值
        t0_transfer = auxSFF.naturalSFF(ii_forma-1).tf;
        tf_transfer = t0_transfer + dt_transfer(ii_forma);
        % 只读最后的点可能会因为浮点数误差读不出来，所以读取两个点，取最后一个点
        x0_j2k_transfer = PosVel_j2k(auxSTK.sat_forma(ii_forma-1).satA.sat,t0_transfer-1,t0_transfer,3600);
        x0_j2k_transfer = x0_j2k_transfer(end,2:end);
    end
    
    % 若satB_sep已存在，则删除原卫星
    % 若transfer卫星不存在，则新建卫星；若存在，则更新初值与递推时间
    satB_trans_name = ['satB_trans',num2str(ii_forma)];
    if  auxSTK.scenario.Children.Contains('eSatellite',satB_trans_name) 
        auxSTK.scenario.Children.Unload('eSatellite',satB_trans_name);
    end
    % 新建transfer卫星
    auxSTK.satB_trans(ii_forma) = sat_initialize(auxSTK.scenario,x0_j2k_transfer,'CentralBody/Earth J2000',...
        t0_transfer,tf_transfer-t0_transfer,satB_trans_name,auxSFF.PropagatorName);
    
    auxSFF.transferSFF(ii_forma).t0 = t0_transfer;
    auxSFF.transferSFF(ii_forma).tf = tf_transfer;
    
    % 重新调整natural卫星B的初始时刻与递推时间
    t0_natural = tf_transfer;
    tf_natural = auxSFF.naturalSFF(ii_forma).tf;
    satA = auxSTK.sat_forma(ii_forma).satA;
    satB = auxSTK.sat_forma(ii_forma).satB;
    
    % 统一读取并存储auxSFF.naturalSFF的数据
    tstep = 3600;
    auxSFF.naturalSFF(ii_forma).t0 = t0_natural;
    auxSFF.naturalSFF(ii_forma) = SaveOrbit2Aux(auxSFF.naturalSFF(ii_forma),...
        satA.sat,satB.sat,t0_natural,tf_natural,tstep);
    
    if ii_forma == 1
        % 存储分离段的轨道数据
        satB_sep = auxSTK.satB_sep;
        sepSFF.t0 = t0_sep;
        sepSFF.tf = tf_sep;
        sepSFF.dv0_j2kLVLH = dv0_sep;
        sepSFF = SaveOrbit2Aux(sepSFF,satA.sat,satB_sep.sat,t0_sep,tf_sep,tstep);
        auxSFF.sepSFF = sepSFF;
    end
end

