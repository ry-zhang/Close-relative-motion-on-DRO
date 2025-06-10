% DE430星历积分

clear; clc;

load('DE430Coeff.mat');
% DE430 = from 14 Dec 1949 00:00:00.000 UTCG (2433264.50000000) 
% to 21 Dec 2149 00:00:00.000 UTCG (2506320.50000000)

%% 全局参数
C_Mat = DE430Coeff;
AU = 0.149597870700000000e+09;
UTC2TCB  = [0,0,0,0, 1, 9.186]; % @2020+之后的跳秒
%赤道系/黄道系转换矩阵
rotation_matrices; 
% 天体引力参数from DE430
xmu(1) = 2.2031780000000021E+04 ; % mercury
xmu(2) = 3.2485859200000006E+05;
xmu(3) = 398600.435436;
xmu(4) = 4.2828375214000022E+04; %mars
xmu(5) = 1.2671276480000021E+08 ; % jupiter
xmu(6) =  3.7940585200000003E+07;   % saturn
xmu(7) = 5.7945486000000080E+06 ;
xmu(8) = 6.8365271005800236E+06 ;
xmu(9) = 9.7700000000000068E+02; %pluto
xmu(10) = 4902.801;  %moon
xmu(11) = 1.3271244004193938E+11; % sun, km^3*s^-2

%% 从UTC计算JD
t0UTC  = [2019 6 14 4 0 0];
t0TCB = t0UTC+UTC2TCB;
jd0 = juliandate(t0TCB(1), t0TCB(2), t0TCB(3), t0TCB(4), t0TCB(5), t0TCB(6));

%% 积分器辅助变量
% 需要配置DE系数/天体引力参数/开始儒略日/坐标转换矩阵/考虑的天体类型
aux.C_Mat = DE430Coeff;         
aux.xmu = xmu;                           
aux.jd0 = jd0; 
aux.ICRF_2_MeanEclpJ2k = ICRF_2_MeanEclpJ2k;
aux.ICRF_2_J2K =ICRF_2_J2K;
aux.threeBodyList = [1:9]; % 考虑的三体类型

%% 例1：太阳系, 在sun Eclp2k 下积分
x0 = [-128392580.265370,-64034354.8862313,-2800775.03434449,...
    15.5056980946092,-22.9424205617426,-0.667381959318738]';
tSpan =  [0 : 2*86400 : 170*86400];

options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, xx_Helio] = ode113(@(t, xx)eqom_helioEclpJ2k(t, xx, aux), tSpan, x0 , options);
try
    clear Satellite1rsun
catch
end

%filename = 'D:\OneDrive\Files\matlab\csu\zhr\helio_EOM\Satellite1_r_sun.txt';
filename = 'Satellite1_r_sun.txt';

startRow = 6;
formatSpec = '%*24s%21f%21f%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
fclose(fileID);
Satellite1rsun = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;

err = xx_Helio(:, 1:3) - Satellite1rsun;
close all;

figure(100) % 积分误差
plot(tSpan/86400, err); hold on; grid on;
xlabel('t [day]')
ylabel('dr [km]')


%% 例2：地月系,  在地心惯性系MEME J2k中积分
x0 = [-26e4 -26e4 -5e4 1.1 -0.8 0]';
tSpan =  [0 : 0.2*86400 : 45*86400];
aux.threeBodyList = [1:2, 4:6]; % 水 金 火 土 木星
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, xx_Geo] = ode113(@(t, xx)eqom_geoMEMEJ2k(t, xx, aux), tSpan, x0 , options);

%filename = 'D:\OneDrive\Files\matlab\csu\zhr\helio_EOM\geo_J2000_Position_Velocity.txt';
filename = 'geo_J2000_Position_Velocity.txt';
startRow = 6;
formatSpec = '%*24s%18f%18f%18f%15f%15f%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);
fclose(fileID);
geoJ2000PositionVelocity = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;

err = xx_Geo - geoJ2000PositionVelocity;

figure(200)
plot(tSpan/86400, err); hold on; grid on;
xlabel('t [day]')
ylabel('dr [km]')

figure(201)
plot3(xx_Geo(:,1), xx_Geo(:,2), xx_Geo(:,3)); hold on; grid on;
axis equal

%%  例3：rv从归一化地月质心旋转系-->地心J2k
rvRot = [0.987849415730060 0 0 0 0 0]';
rvJ2k  = BCRot_To_ECJ2k(jd0, rvRot, DE430Coeff);
rvJ2k = [-248516, -226961, -67695, 0.981004, -0.93483, -0.467167]';
rvRot1  = ECJ2k_To_ECRot(jd0, rvJ2k, DE430Coeff);
rvRot2  = ECJ2k_To_MCRot(jd0, rvJ2k, DE430Coeff);

%% 例4：地月DRO积分例子，原始数据@归一化地月质心旋转系
%filename = 'D:\OneDrive\Files\matlab\csu\zhr\helio_EOM\DRO.txt';
filename = 'DRO.txt';

formatSpec = '%22f%23f%24f%23f%23f%23f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
DRO = [dataArray{1:end-1}];
clearvars filename formatSpec fileID dataArray ans;

for ii = 1:3000:size(DRO,1) % 挑几条算一算
    rvRot =  DRO(ii,1:6)';    
    x0  = BCRot_To_ECJ2k(jd0, rvRot, DE430Coeff);

    tSpan =  [0 : 0.2*86400 : 100*86400];
    aux.threeBodyList = [];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~, xx_Dro] = ode113(@(t, xx)eqom_geoMEMEJ2k(t, xx, aux), tSpan, x0 , options);   
    figure(402)
    plot3(xx_Dro(:,1), xx_Dro(:,2), xx_Dro(:,3)); hold on; grid on;
    axis equal
end

%% 例5：地月周期轨道积分后在不同坐标系展示，原始数据@归一化地月质心旋转系
rvRot =  DRO(9001,1:6)'; % DRO初值
rvRot =  [0.824168945125376, 0, 0.022302875780033, 0, 0.133751712663525, 0]'; % NRHO初值
rvRot = [1.11792521878049,0,0.0183504776974601,0,0.183156804395910,0]';  %Halo初值
x0  = BCRot_To_ECJ2k(jd0, rvRot, DE430Coeff);
tSpan =  [0 : 0.1*86400 : 7*86400];
aux.threeBodyList = [];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, xx_ECI] = ode113(@(t, xx)eqom_geoMEMEJ2k(t, xx, aux), tSpan, x0 , options);

xx_ECRot = [];  xx_MCRot = [];
for ii = 1:length(t)
    xx_ECRot(ii, :) = (ECJ2k_To_ECRot(jd0+t(ii)/86400, xx_ECI(ii, :)', DE430Coeff))';
    xx_MCRot(ii, :) = (ECJ2k_To_MCRot(jd0+t(ii)/86400, xx_ECI(ii, :)', DE430Coeff))';
end

close all
figure(501)
subplot(2,2,1)
plot3(xx_ECI(:,1), xx_ECI(:,2), xx_ECI(:,3)); hold on; grid on;
axis equal
xlabel('x');  ylabel('y');  zlabel('z'); 
subplot(2,2,2)
plot3(xx_ECRot(:,1), xx_ECRot(:,2), xx_ECRot(:,3)); hold on; grid on;
axis equal
xlabel('x');  ylabel('y');  zlabel('z'); 
subplot(2,2,3)
plot3(xx_MCRot(:,1), xx_MCRot(:,2), xx_MCRot(:,3)); hold on; grid on;
plot3(xx_MCRot(1,1), xx_MCRot(1,2), xx_MCRot(1,3), 'ro'); hold on; grid on;
axis equal
xlabel('x');  ylabel('y');  zlabel('z'); 