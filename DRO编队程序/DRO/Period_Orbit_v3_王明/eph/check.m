
function check
clc; clear; close all; clear global;
addpath('D:\360安全云盘同步版\[Code]\func');
% addpath('D:\360安全云盘同步版\[Code]\jpl_ephem');

% addpath('de421');
% addpath('de421'); % 使用de421星历
% global iephem ephname km
% iephem = 1;
% ephname = 'de421.bin';
% km = 1;

load DE430_1949to2150
load DE430

dd = [];

% save DE430_1949to2150 DE430_1949to2150
jd = jday(2020 , 1 , 1 , 0 , 0 , 0);
ntarg = 3;
ncent = 11;

% load Coeff2020to2040
% Coeff2020to2040 = DE430Coeff;
% Coeff2020to2040 = Coeff2020to2040(800 : 1028 , :);
% save Coeff2020to2040 Coeff2020to2040

% tic
% for iLoop = 1 : 1e3
%     rrd = jplephem (jd, ntarg, ncent);
% end
% toc

% [mon,day,year] = gdate(DE430Coeff(1 , 1))
% [mon,day,year] = gdate(DE430Coeff(end , 2))

tic
for iLoop = 1 : 1e5
    
%     rv = ephEme(jd , ntarg, ncent , DE430); % 3.5
%         rv = ephEme_mex(jd , ntarg, ncent , DE430); % 2.0
    
%     rv = ephEclip(jd , ntarg, ncent , DE430); % 3.6
%         rv = ephEclip_mex(jd , ntarg, ncent , DE430); % 2.3
    
end
toc

dd = [];

end
