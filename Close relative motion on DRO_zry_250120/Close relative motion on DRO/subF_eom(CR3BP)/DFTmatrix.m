function D = DFTmatrix(N)
% matrix of FFT
%
% D = DFTmatrix(N)
%
% Input arguments:
% -------------------------------------------------------------------------
% N          [1x1]               FFT size, interger
% 
% Output arguments:
% -------------------------------------------------------------------------
% D          [NxN]               Matrix of FFT
% 
% External functions called:
% -------------------------------------------------------------------------
%  none
% 
% Copyright (C) 25/7/2020 by Chihang Yang 
% email: ychhtl@foxmail.com

%%
k_all = linspace(floor(-(N-1)/2),floor((N-1)/2),N)';
% k_all = linspace(0,N-1,N)';
j_all = linspace(0,N-1,N);
D = exp(-1i*k_all*2*pi*j_all/N);