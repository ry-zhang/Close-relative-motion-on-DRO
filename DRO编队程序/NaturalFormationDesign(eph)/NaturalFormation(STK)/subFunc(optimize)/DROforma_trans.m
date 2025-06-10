function [f,dv2] = DROforma_trans(dv0_scaled,x0_j2kLVLH_REL,xf_j2kLVLH_REL_tar,scale_factor,satB,satA)
% 脉冲变轨约束函数

dv0 = dv0_scaled./scale_factor;
x0_j2kLVLH_REL(4:6) = x0_j2kLVLH_REL(4:6) + dv0;

% 赋值给satB的初值
satA_Name = satA.sat.InstanceName;
satB.InitialState.CoordSystemName = ['Satellite/',satA_Name,' J2KLVLH_SC'];
satB.InitialState.SetElementType('eVAElementTypeCartesian');
set(satB.InitialState.Element, 'X', x0_j2kLVLH_REL(1),'Y',x0_j2kLVLH_REL(2),...
    'Z',x0_j2kLVLH_REL(3),'Vx',x0_j2kLVLH_REL(4),...
    'Vy',x0_j2kLVLH_REL(5),'Vz',x0_j2kLVLH_REL(6));
satB.sat.Propagator.RunMCS;

% 积分，计算相对运动rv_TC_LVLH
t0 = satB.InitialState.OrbitEpoch;
tf = t0 + satB.Propagate.StoppingConditions.Item('Duration').Properties.Trip;
tstep = 3600;
rf_j2kLVLH_REL = PosVel_j2kLVLH(satB.sat,satA.sat,tf-1,tf,tstep);
rf_j2kLVLH_REL = rf_j2kLVLH_REL(end,:);

% 末端位置误差
f = rf_j2kLVLH_REL(1:3) - xf_j2kLVLH_REL_tar(1:3);

% 末端速度误差
dv2 = xf_j2kLVLH_REL_tar(4:6) - rf_j2kLVLH_REL(4:6);