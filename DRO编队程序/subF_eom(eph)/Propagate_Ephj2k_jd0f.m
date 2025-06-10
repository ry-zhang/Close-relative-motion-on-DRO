function xf_j2k = Propagate_Ephj2k_jd0f(x0_j2k,jd0,jdf,DynamicModel)
% 初值
% jd0 = 2460109.81760632; % TBD
% x0_j2k = [260861.357166097,167404.871364748,78987.9106263715,...
%         -0.797706314601184,0.975263733684134,0.530154969412835];
% jdf = juliandate(datetime('2023-01-01') + minutes(1) + seconds(9.186)); % TBD

% 将初值积分到期望时刻
aux = []; 
load('DE430Coeff.mat');%星历表
aux.C_Mat = DE430Coeff;
aux.jd0 = jd0; % TBD

if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,1); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux = initialize(aux); % 初始化
else
    error('Wrong DynamicModel')
end

tspan_sec = [0,(jdf-aux.jd0)*86400]; % TBD
xf_j2k = Propagate_Ephj2k(x0_j2k,tspan_sec,tspan_sec(2),aux);