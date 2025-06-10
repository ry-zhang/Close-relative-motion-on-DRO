function [auxSFF,auxSTK] = sat_naturalSFF(auxSFF,auxSTK,dt_all_0)
% 计算入轨之后抵达theta_tar的滑行时间，在STK中建立编队卫星satA_form与satB_form
T_Moon = (23*60+56-52.6)*60; %月球上中天的大致周期
for ii_forma = 1:length(dt_all_0)
    % 设置forma段的初始数据
    if ii_forma == 1
        auxSFF.naturalSFF(ii_forma).t0 = auxSFF.coasting.tf;
        auxSFF.naturalSFF(ii_forma).x0_MCEMR_A = auxSFF.coasting.xx_MCEMR_AB(end,:);
        auxSFF.naturalSFF(ii_forma).x0_j2k_A = auxSFF.coasting.xx_j2k_AB(end,:);
    else
        auxSFF.naturalSFF(ii_forma).t0 = auxSFF.naturalSFF(ii_forma-1).tf;
        auxSFF.naturalSFF(ii_forma).x0_MCEMR_A = auxSFF.naturalSFF(ii_forma-1).xx_MCEMR_A(end,:);
        auxSFF.naturalSFF(ii_forma).x0_j2k_A = auxSFF.naturalSFF(ii_forma-1).xx_j2k_A(end,:);
    end
    
    % 创建A星
    % 每一个卫星创建两段propagate，第一段递推时间为dt，第二段是20day，然后计算第二段至theta_tar的具体时间
    x0_j2k = auxSFF.naturalSFF(ii_forma).x0_j2k_A;
    t0_iniState = auxSFF.naturalSFF(ii_forma).t0;
    dt_p1 = dt_all_0(ii_forma);
    dt_p2_0 = 20*T_Moon; % 明显大于DRO的周期即可，确保星历下足够绕月一周即可
    sat_A_form = sat_initialize_forma(auxSTK.scenario,...
        x0_j2k,t0_iniState,dt_p1,dt_p2_0,['sat_A_form',num2str(ii_forma)],auxSFF.PropagatorName);
    auxSTK.sat_forma(ii_forma).satA = sat_A_form;
    
    % 计算A星的轨道
    t0_form1_sec = sat_A_form.InitialState.OrbitEpoch;
    tf1_form1_sec = t0_form1_sec + sat_A_form.Propagate1.StoppingConditions.Item('Duration').Properties.Trip;
    tf2_form1_sec_0 = tf1_form1_sec + sat_A_form.Propagate2.StoppingConditions.Item('Duration').Properties.Trip;
    
%     xx_A_form1_MCEMR = PosVel_MCEMR(sat_A_form.sat,tf1_form1_sec,tf2_form1_sec_0,3600);
    % 将间隔调成月球上中天的时间,以保证每次变轨的时间在地球可观测区内
    xx_A_form1_MCEMR = PosVel_MCEMR(sat_A_form.sat,tf1_form1_sec,tf2_form1_sec_0,T_Moon);
    % 判断在tf1_form1_sec至tf2_form1_sec_0之间，抵达theta_tarAll(1)与theta_tarAll(2)的时间，选取其中较小的，作为目标时间
    [tf_form1_sec1,xf_form1_MCEMR1,tf_form1_index1] = find_thetatar(xx_A_form1_MCEMR,auxSFF.theta_tarAll(1));
    [tf_form1_sec2,xf_form1_MCEMR2,tf_form1_index2] = find_thetatar(xx_A_form1_MCEMR,auxSFF.theta_tarAll(2));
    if tf_form1_sec1<tf_form1_sec2
        auxSFF.alpha_tar(ii_forma+1) = auxSFF.alpha_tarAll(1);
        auxSFF.theta_tar(ii_forma+1) = auxSFF.theta_tarAll(1);
        tf_form1_sec = tf_form1_sec1;
    else
        auxSFF.alpha_tar(ii_forma+1) = auxSFF.alpha_tarAll(2);
        auxSFF.theta_tar(ii_forma+1) = auxSFF.theta_tarAll(2);
        tf_form1_sec = tf_form1_sec2;
    end
    sat_A_form.Propagate2.StoppingConditions.Item('Duration').Properties.Trip = tf_form1_sec - tf1_form1_sec;
    sat_A_form.sat.Propagator.RunMCS;
    xx_A_form1_MCEMR = PosVel_MCEMR(sat_A_form.sat,t0_form1_sec,tf_form1_sec,3600);
    xx_A_form1_j2k = PosVel_j2k(sat_A_form.sat,t0_form1_sec,tf_form1_sec,3600);

    % 保存format段的末端数据
    auxSFF.naturalSFF(ii_forma).tf = tf_form1_sec;
    auxSFF.naturalSFF(ii_forma).t_sample = xx_A_form1_MCEMR(:,1);
    auxSFF.naturalSFF(ii_forma).xx_MCEMR_A = xx_A_form1_MCEMR(:,2:end);
    auxSFF.naturalSFF(ii_forma).xx_j2k_A = xx_A_form1_j2k(:,2:end);
    
    % 按照A星的数据建造B星，以便后续优化
    dt_p2 = tf_form1_sec - tf1_form1_sec;
    sat_B_form = sat_initialize_forma(auxSTK.scenario,...
        x0_j2k,t0_iniState,dt_p1,dt_p2,['sat_B_form',num2str(ii_forma)],auxSFF.PropagatorName);
    auxSTK.sat_forma(ii_forma).satB = sat_B_form;
end

