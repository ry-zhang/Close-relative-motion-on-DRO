function PosVel = PosVel_j2k(sat,t0,tf,tstep)
% 读取MCEMR坐标系中的航天器的位置速度

rptElemsPos = {'Time';'x';'y';'z'}; 
repPos = sat.DataProviders.Item('Cartesian Position');
repPos_aa = repPos.Group.Item('J2000');
DataIntePos = repPos_aa.ExecElements(t0,tf,tstep,rptElemsPos);
Pos = cell2mat(DataIntePos.DataSets.ToArray);

rptElemsVel = {'x';'y';'z'}; 
repVel = sat.DataProviders.Item('Cartesian Velocity');
repVel_aa = repVel.Group.Item('J2000');
DataInteVel = repVel_aa.ExecElements(t0,tf,tstep,rptElemsVel);

Vel = cell2mat(DataInteVel.DataSets.ToArray);
PosVel = [Pos,Vel];

