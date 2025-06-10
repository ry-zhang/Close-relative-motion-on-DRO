function rvMC_MCR_c = T_LVLH2MCR(rvTC_LVLH_c,rvMT_MCR_t)
%�������ϵLVLH���ǹ�һ����ת�������ģ�����ԭ�㣩��תϵMCR(��һ����
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH
% LVLH����ϵ��RIC

%���룺
%    Ŀ�꺽����λ���ٶ�rMT_MCR,vMT_MCR
%    ���ٺ�����λ���ٶ�rTC_LVLH,vTC_LVLH
%�����
%    ���ٺ�����λ���ٶ�rTC_MCR,vTC_MCR

%�������ϵLVLH
%    -ԭ��λ��Ŀ�꺽��������
%    -x��������������ָ��Ŀ�꺽����
%    -z������Ŀ�꺽�����������������Ľ��ٶȷ��� 
%������תϵMCR
%    -ԭ��λ������
%    -x��������ָ�������
%    -z��Ϊ˲ʱ��������ڵ���ĽǶ�������
row = size(rvTC_LVLH_c,1);
rvMC_MCR_c = zeros(row,6);
for ii = 1:row
    rTC_LVLH_c_ii = rvTC_LVLH_c(ii,1:3); vTC_LVLH_c_ii = rvTC_LVLH_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_t_ii = rvMT_MCR_t(ii,4:6);
    [rTC_MCR_c_ii,vTC_MCR_c_ii] = LVLH2MCR(rTC_LVLH_c_ii',vTC_LVLH_c_ii',rMT_MCR_t_ii',vMT_MCR_t_ii');
    rvMC_MCR_c(ii,:) = [rTC_MCR_c_ii',vTC_MCR_c_ii'];
end

end

function [rMC_MCR_c,vMC_MCR_c] = LVLH2MCR(rTC_LVLH_c,vTC_LVLH_c,rMT_MCR_t,vMT_MCR_t)
hT = cross(rMT_MCR_t,vMT_MCR_t);
hT_norm = norm(hT);
rMT_norm = norm(rMT_MCR_t);
omegaT = hT/rMT_norm^2;
%LVLH����ϵ����
iLVLH = rMT_MCR_t/rMT_norm;
kLVLH = hT/hT_norm;
jLVLH = cross(kLVLH,iLVLH);

%����ת��
TM_LVLH2MCR = [iLVLH jLVLH kLVLH];
rTC_MCR = TM_LVLH2MCR*rTC_LVLH_c;
rMC_MCR_c = rTC_MCR;
% rMC_MCR_c = rMT_MCR_t+rTC_MCR;

vTC_MCR = TM_LVLH2MCR*vTC_LVLH_c+cross(omegaT,rTC_MCR);
vMC_MCR_c = vTC_MCR;
% vMC_MCR_c = vMT_MCR_t+vTC_MCR;
end

