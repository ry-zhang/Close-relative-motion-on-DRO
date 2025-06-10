function PosVel = PosVel_MCEMR(sat,t0,tf,tstep)
% 读取MCEMR坐标系中的航天器的位置速度

rptElems = {'Time';'x';'y';'z';'Velocity x';'Velocity y';'Velocity z'}; 
repPos = sat.DataProviders.Item('Points Choose System');
repPos_aa = repPos.Group.Item('Center');
repPos_aa.PreData = 'CentralBody/Moon MCEMR';
DataInte = repPos_aa.ExecElements(t0,tf,tstep,rptElems);

PosVel = cell2mat(DataInte.DataSets.ToArray);

