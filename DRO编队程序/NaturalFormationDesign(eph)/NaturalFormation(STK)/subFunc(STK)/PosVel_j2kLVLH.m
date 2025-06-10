function PosVel = PosVel_j2kLVLH(satB,satA,t0,tf,tstep)
% 读取MCEMRLVLH坐标系中的航天器的位置速度
satA_Name = satA.InstanceName;

rptElems = {'x';'y';'z';'Velocity x';'Velocity y';'Velocity z'}; 
repPos = satB.DataProviders.Item('Points Choose System');
repPos_aa = repPos.Group.Item('Center');
repPos_aa.PreData = ['Satellite/',satA_Name,' J2KLVLH_SC'];
DataInte = repPos_aa.ExecElements(t0,tf,tstep,rptElems);

PosVel = cell2mat(DataInte.DataSets.ToArray);

