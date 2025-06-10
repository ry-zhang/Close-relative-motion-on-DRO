function f = DROforma_dist(x0,rr_MCEMRLVLH_REL_ref,scale_factor,satB,satA)
% 任务轨道优化目标函数

x0 = x0./scale_factor;
x0_MCEMRLVLH_REL = [x0(1),x0(2),0,x0(3),x0(4),0];

% 赋值给satB的初值
satA_Name = satA.sat.InstanceName;
satB.InitialState.CoordSystemName = ['Satellite/',satA_Name,' MCEMRLVLH_SC'];
satB.InitialState.SetElementType('eVAElementTypeCartesian');
set(satB.InitialState.Element, 'X', x0_MCEMRLVLH_REL(1),'Y',x0_MCEMRLVLH_REL(2),...
    'Z',x0_MCEMRLVLH_REL(3),'Vx',x0_MCEMRLVLH_REL(4),...
    'Vy',x0_MCEMRLVLH_REL(5),'Vz',x0_MCEMRLVLH_REL(6));
satB.sat.Propagator.RunMCS;

% 积分，计算相对运动rv_TC_LVLH
t0 = satB.InitialState.OrbitEpoch;
tf_p1 = t0 + satB.Propagate1.StoppingConditions.Item('Duration').Properties.Trip;
tf = tf_p1 + satB.Propagate2.StoppingConditions.Item('Duration').Properties.Trip;
tstep = 3600;
rr_MCEMRLVLH_REL = PosVel_MCEMRLVLH(satB.sat,satA.sat,t0,tf,tstep);

r_error = sqrt(sum((rr_MCEMRLVLH_REL(:,[1,2,3])-rr_MCEMRLVLH_REL_ref).^2,2));
% 距离参考轨道的最大距离
f = max(r_error);
% r_norm = sqrt(sum((rr_MCEMRLVLH_REL_ref).^2,2));
% f = max(r_error./r_norm*max(r_norm));