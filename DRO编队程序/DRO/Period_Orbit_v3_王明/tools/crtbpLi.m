function [Li_pos , Li_c] = crtbpLi(mu)
% ����crtbp��5���������յ������
%
% ���ӣ�
% % [1]
% aux = crtbpAuxEM;
% [Li_pos , Li_c] = crtbpLi(aux.mu)
%
% % [2]
% aux = crtbpAuxSE;
% [Li_pos , Li_c] = crtbpLi(aux.mu)
%
% 2015/4/10
% Copyright(C) Chen Zhang
% -------------------------------------------------------------------------
% ���ó�ʼ�²�ֵ
x012 = (mu/3)^(1/3) ;
x03  = 1-(7/12)*mu  ;
% ����5��ƽ��������
options = optimset('TolFun', 2.5e-14, 'TolX', 2.5e-14) ;
L1  = 1 - mu - fzero(@(x) x^5-(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2+2*mu*x-mu, ...
    x012, options) ;
L2  = 1 - mu + fzero(@(x) x^5+(3-mu)*x^4+(3-2*mu)*x^3-mu*x^2-2*mu*x-mu, ...
    x012, options) ;
L3  =   - mu - fzero(@(x) x^5+(2+mu)*x^4+(1+2*mu)*x^3-(1-mu)*x^2-2*(1-mu)*x...
    -(1-mu), x03,  options) ;
L4x = 0.5 - mu ;
L4y = sqrt(3)/2 ;
% ����5��ƽ����λ��
Li_pos = [L1, 0;
    L2, 0;
    L3, 0;
    L4x, L4y;
    L4x, -L4y];
% ����������1�����5��ƽ��������
if nargout > 1
    xx = [Li_pos , zeros(5 , 2)];
    Li_c = crtbpJacobi2D(xx , mu);
end
end
