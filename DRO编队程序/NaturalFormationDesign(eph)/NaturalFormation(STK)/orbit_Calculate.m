% 根据优化结果，重新创建STK工程文件，计算轨道数据
load('DRO_SSF.mat')

%% Launch a new instance of STK11 and grab it
stkMode = 1;
if stkMode ==1
    % 方法1，可以控制stk程序是否可见,调试方便
    uiapp = actxserver('STK11.application');
    uiapp.Visible = 1;
    root = uiapp.Personality2;
elseif stkMode ==2
    % 方法2，stk在后台运行,速度快
    STKXApplication = actxserver('STKX11.application');
    root = actxserver('AgStkObjects11.AgStkObjectRoot');
end

% create a new scenario
scenario = root.Children.New('eScenario','DROFormation');
% set the time period
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('JED'); 
scenario.SetTimePeriod(num2str(auxSFF.jd0), ['+', num2str(365), ' days']);
% reset the animation
root.ExecuteCommand('Animate * Reset');
% set the unit of time as Epoch second to convinient the calculation
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec'); 
% 创建月心旋转坐标系
CreateMCRframe(scenario,root)
PropagatorName = 'Earth HPOP Default v10';

%% 建立AB两颗卫星
% ------------------------------建立A星------------------------------------
satA = sat_initialize(scenario,auxSFF.x00_j2k_AB, 'CentralBody/Earth J2000',...
    scenario.StartTime, auxSFF.naturalSFF(end).tf, 'satA',PropagatorName);
auxSAT.satA = satA;
% 建立A星相对运动坐标系
CreateMCRLVLHframe(auxSAT);
CreateJ2KLVLHframe(auxSAT);

% ------------------------------建立B星------------------------------------
% 初始状态，滑行段，分离脉冲，分离段
% dv0_1，变轨段1，dvf_1,任务段1
% dv0_2，变轨段2，dvf_2,任务段2
% dv0_3，变轨段3，dvf_3,任务段3

satB.sat = scenario.Children.New('eSatellite','satB');
satB.sat.SetPropagatorType('ePropagatorAstrogator');
satBMSC = satB.sat.Propagator.MainSequence;
satBMSC.RemoveAll;

% insert and set the initial state
satB.InitialState = satBMSC.Insert('eVASegmentTypeInitialstate','InitialState','-');
satB.InitialState.OrbitEpoch = scenario.StartTime;
satB.InitialState.CoordSystemName = 'CentralBody/Earth J2000';
satB.InitialState.SetElementType('eVAElementTypeCartesian');
set(satB.InitialState.Element, 'X', auxSFF.x00_j2k_AB(1),'Y',auxSFF.x00_j2k_AB(2),...
    'Z',auxSFF.x00_j2k_AB(3),'Vx',auxSFF.x00_j2k_AB(4),'Vy',auxSFF.x00_j2k_AB(5),'Vz',auxSFF.x00_j2k_AB(6));
% insert and set the Coast segment
satB.Coast = satBMSC.Insert('eVASegmentTypePropagate','Coast','-');
satB.Coast.PropagatorName = PropagatorName;
satB.Coast.StoppingConditions.Item('Duration').Properties.Trip = auxSFF.coasting.tf-auxSFF.coasting.t0;
% insert and set the seperation Maneuver
satB.SepManeuver = satBMSC.Insert('eVASegmentTypeManeuver','SepManeuver','-');
satB.SepManeuver.SetManeuverType('eVAManeuverTypeImpulsive');
satB.SepManeuver.Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
satB.SepManeuver.Maneuver.AttitudeControl.ThrustAxesName = 'Satellite/satA J2KLVLH';
dv0_sep = auxSFF.sepSFF.dv0_j2kLVLH;
satB.SepManeuver.Maneuver.AttitudeControl.DeltaVVector.AssignCartesian(dv0_sep(1),dv0_sep(2),dv0_sep(3));
% insert and set the seperation Propagate
satB.SepPropagate = satBMSC.Insert('eVASegmentTypePropagate','SepPropagate','-');
satB.SepPropagate.PropagatorName = PropagatorName;
satB.SepPropagate.StoppingConditions.Item('Duration').Properties.Trip = auxSFF.sepSFF.tf-auxSFF.sepSFF.t0;

for ii_forma = 1:length(auxSFF.naturalSFF)
    % insert and set the formation dv0
    satB.Maneuver0(ii_forma) = satBMSC.Insert('eVASegmentTypeManeuver',['Maneuver0',num2str(ii_forma)],'-');
    satB.Maneuver0(ii_forma).SetManeuverType('eVAManeuverTypeImpulsive');
    satB.Maneuver0(ii_forma).Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
    satB.Maneuver0(ii_forma).Maneuver.AttitudeControl.ThrustAxesName = 'Satellite/satA J2KLVLH';
    dv0_ii = auxSFF.transferSFF(ii_forma).dv0_j2kLVLH;
    satB.Maneuver0(ii_forma).Maneuver.AttitudeControl.DeltaVVector.AssignCartesian(dv0_ii(1),dv0_ii(2),dv0_ii(3));
    
    % insert and set the formation Maneuver propagate segment
    satB.Propagate_transf(ii_forma) = satBMSC.Insert('eVASegmentTypePropagate',['Propagate_transf',num2str(ii_forma)],'-');
    satB.Propagate_transf(ii_forma).PropagatorName = PropagatorName;
    satB.Propagate_transf(ii_forma).StoppingConditions.Item('Duration').Properties.Trip = auxSFF.transferSFF(ii_forma).tf-auxSFF.transferSFF(ii_forma).t0;
    
    % insert and set the formation dvf
    satB.Maneuverf(ii_forma) = satBMSC.Insert('eVASegmentTypeManeuver',['Maneuverf',num2str(ii_forma)],'-');
    satB.Maneuverf(ii_forma).SetManeuverType('eVAManeuverTypeImpulsive');
    satB.Maneuverf(ii_forma).Maneuver.SetAttitudeControlType('eVAAttitudeControlThrustVector');
    satB.Maneuverf(ii_forma).Maneuver.AttitudeControl.ThrustAxesName = 'Satellite/satA J2KLVLH';
    dvf_ii = auxSFF.transferSFF(ii_forma).dvf_j2kLVLH;
    satB.Maneuverf(ii_forma).Maneuver.AttitudeControl.DeltaVVector.AssignCartesian(dvf_ii(1),dvf_ii(2),dvf_ii(3));
    
    % insert and set the formation Maneuver0 propagate segment
    satB.Propagate(ii_forma) = satBMSC.Insert('eVASegmentTypePropagate',['Propagate',num2str(ii_forma)],'-');
    satB.Propagate(ii_forma).PropagatorName = PropagatorName;
    satB.Propagate(ii_forma).StoppingConditions.Item('Duration').Properties.Trip = auxSFF.naturalSFF(ii_forma).tf-auxSFF.naturalSFF(ii_forma).t0;
end
% run the astrogator
satB.sat.Propagator.RunMCS;

%% 按段读取轨道并保存数据
tstep = 60; % 采样点间隔，sec
auxSFF.sepSFF = SaveOrbit2Aux(auxSFF.sepSFF,satA.sat,satB.sat,auxSFF.sepSFF.t0,auxSFF.sepSFF.tf,tstep);
auxSFF.coasting = SaveOrbit2Aux(auxSFF.coasting,satA.sat,satB.sat,auxSFF.coasting.t0,auxSFF.coasting.tf,tstep);
for ii_forma = 1:length(auxSFF.naturalSFF)
    auxSFF.naturalSFF(ii_forma) = SaveOrbit2Aux(auxSFF.naturalSFF(ii_forma),...
        satA.sat,satB.sat,auxSFF.naturalSFF(ii_forma).t0,auxSFF.naturalSFF(ii_forma).tf,tstep);
    auxSFF.transferSFF(ii_forma) = SaveOrbit2Aux(auxSFF.transferSFF(ii_forma),...
        satA.sat,satB.sat,auxSFF.transferSFF(ii_forma).t0,auxSFF.transferSFF(ii_forma).tf,tstep);
end

% 保存数据
save('DRO_SSF_60sec','auxSFF')
