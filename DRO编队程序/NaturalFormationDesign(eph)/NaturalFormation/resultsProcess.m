dbstop if error

clear
% close all
addpath('../../subF_eom(CR3BP)')
addpath('../../subF_eom(eph)')

format longg; format compact
warning off

set(0,'defaultAxesFontName', 'TimesSimSun');%坐标轴
set(0,'defaultTextFontName', 'TimesSimSun');%文字

%%
load('PhDNaturalSSF_DR03_HighFidelity.mat')
naturalSFF = auxSFF.naturalSFF;
ii_segment = 1;

aux = []; 
load('DE430Coeff.mat'); aux.C_Mat = DE430Coeff;
aux.jd0 = auxSFF.jd0 + naturalSFF(ii_segment).t0/86400;
if strcmp(DynamicModel,'HighFidelity')
    aux = SPICEinitialize(aux,1); % 初始化
elseif strcmp(DynamicModel,'SEMephemeris')
    aux = initialize(aux); % 初始化
else
    error('Wrong DynamicModel')
end
x0_j2k_chief = naturalSFF(ii_segment).x0_j2k_target;
x0_j2kLVLH_rel = naturalSFF(ii_segment).x0_j2kLVLH_rel;
a0_j2k_chief = eomj2kMtx(x0_j2k_chief,0,aux);
x0_j2k_deputy = T_TCO2TCR_eph(x0_j2kLVLH_rel,x0_j2k_chief,a0_j2k_chief,'LVLH')+x0_j2k_chief;
% t0UTC = aux.t0UTC;
% jd0 = juliandate(aux.t0UTC);
% datetime(jd0,'ConvertFrom','juliandate','Format','yyyy-MM-dd HH:mm:ss.SSSSSS');
aux = rmfield(aux,'C_Mat');
% save DROformation aux x0_j2k_chief x0_j2k_deputy