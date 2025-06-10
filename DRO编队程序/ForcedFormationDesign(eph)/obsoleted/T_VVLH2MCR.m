function rvMC_MCR_c = T_VVLH2MCR(rvTC_VVLH_c,rvMT_MCR_t)
%�������ϵVVLH���ǹ�һ����ת�������ģ�����ԭ�㣩��תϵMCR(��һ����
% v1 2020/8/20 fhl
% v2 2021/5/28 YCH
%���룺
%Ŀ�꺽����λ���ٶ�rMT_MCR,vMT_MCR
%���ٺ�����λ���ٶ�rTC_VVLH,vTC_VVLH
%�����
%���ٺ�����λ���ٶ�rTC_MCR,vTC_MCR
%�������ϵVVLH
%    -ԭ��λ��Ŀ�꺽��������
%    -y������Ŀ�꺽�����������������Ľ��ٶȵĸ�����
%    -z������Ŀ�꺽����ָ����������
%������תϵMCR
%    -ԭ��λ������
%    -x��������ָ�������
%    -z��Ϊ˲ʱ��������ڵ���ĽǶ�������
row = size(rvTC_VVLH_c,1);
rvMC_MCR_c = zeros(row,6);
for ii = 1:row
    rTC_VVLH_c_ii = rvTC_VVLH_c(ii,1:3); vTC_VVLH_c_ii = rvTC_VVLH_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_t_ii = rvMT_MCR_t(ii,4:6);
    [rTC_MCR_c_ii,vTC_MCR_c_ii] = VVLH2MCR(rTC_VVLH_c_ii',vTC_VVLH_c_ii',rMT_MCR_t_ii',vMT_MCR_t_ii');
    rvMC_MCR_c(ii,:) = [rTC_MCR_c_ii',vTC_MCR_c_ii'];
end

end

function [rMC_MCR_c,vMC_MCR_c] = VVLH2MCR(rTC_VVLH_c,vTC_VVLH_c,rMT_MCR_t,vMT_MCR_t)
hT = cross(rMT_MCR_t,vMT_MCR_t);
hT_norm = norm(hT);
rMT_norm = norm(rMT_MCR_t);
omegaT = hT/rMT_norm^2;
%VVLH����ϵ����
kVVLH = -rMT_MCR_t/rMT_norm;
jVVLH = -hT/hT_norm;
iVVLH = cross(jVVLH,kVVLH);

%����ת��
TM_VVLH2MCR = [iVVLH jVVLH kVVLH];
rTC_MCR = TM_VVLH2MCR*rTC_VVLH_c;
rMC_MCR_c = rTC_MCR;
% rMC_MCR_c = rMT_MCR_t+rTC_MCR;

vTC_MCR = TM_VVLH2MCR*vTC_VVLH_c+cross(omegaT,rTC_MCR);
vMC_MCR_c = vTC_MCR;
% vMC_MCR_c = vMT_MCR_t+vTC_MCR;
end

