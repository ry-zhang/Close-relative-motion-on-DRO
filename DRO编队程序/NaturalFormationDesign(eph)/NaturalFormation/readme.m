% 数据结构说明

% aux_DynModel:         动力学模型参数

% aux_SFF:                  航天器轨道数据
% aux_SFF.jd0UTC            初始时刻的儒略日(UTC)
% aux_SFF.x00_j2k_A    初始时刻的A星状态
% 
% aux_SFF.naturalSFF        三段自然绕飞轨道的始末时刻(为从jd0出发的秒数) 与 A星、B星、二者相对运动的轨道数据(B星-A星)(km,km/s)
% aux_SFF.naturalSFF(1)     短距离伴飞段
% aux_SFF.naturalSFF(2)     中等距离伴飞段
% aux_SFF.naturalSFF(3)     远距离伴飞段
% 
% aux_SFF.transferSFF       四段变轨轨道的始末时刻(为从jd0出发的秒数) 与 始末变轨脉冲(m/s)
% aux_SFF.transferSFF(1)    分离段
% aux_SFF.transferSFF(2)    分离段末端 至 短距离伴飞段的转移
% aux_SFF.transferSFF(3)    短距离伴飞段 至 中等距离伴飞段的转移
% aux_SFF.transferSFF(4)    中等距离伴飞段 至 远距离伴飞段的转移