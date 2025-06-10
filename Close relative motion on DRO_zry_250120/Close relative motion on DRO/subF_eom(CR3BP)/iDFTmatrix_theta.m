function iD = iDFTmatrix_theta(N,theta_all)
% matrix of inverse FFT with frequency shift -rho
%
% iD = iDFTmatrix(N,rho)
%
% Input arguments:
% -------------------------------------------------------------------------
% N          [1x1]               FFT size, interger
% theta_all  [1xN2]              target theta
% 
% Output arguments:
% -------------------------------------------------------------------------
% iD         [NxN]               Matrix of inverse FFT
% 
% External functions called:
% -------------------------------------------------------------------------
%  none
% 
% Copyright (C) 25/7/2020 by Chihang Yang 
% email: ychhtl@foxmail.com

%%
theta_all = theta_all(:);
k_all = linspace(floor(-(N-1)/2),floor((N-1)/2),N);
% k_all = linspace(0,N-1,N);
iD = exp(1i*(theta_all)*k_all)/N;