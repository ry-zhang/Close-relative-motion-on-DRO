function rvMC_MCR_c = T_VNC2MCR(rvTC_VNC_c,rvMT_MCR_t)
%�������ϵVNC���ǹ�һ����ת�������ģ�����ԭ�㣩��תϵMCR(��һ����
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH

%���룺
%    Ŀ�꺽����λ���ٶ�rMT_MCR,vMT_MCR
%    ���ٺ�����λ���ٶ�rTC_VNC,vTC_VNC
%�����
%    ���ٺ�����λ���ٶ�rTC_MCR,vTC_MCR

%�������ϵVNC
%    -ԭ��λ��Ŀ�꺽��������
%    -x������Ŀ�꺽�������ٶȷ���
%    -y������Ŀ�꺽�����������������Ľ��ٶȷ��� 
%������תϵMCR
%    -ԭ��λ������
%    -x��������ָ�������
%    -z��Ϊ˲ʱ��������ڵ���ĽǶ�������
row = size(rvTC_VNC_c,1);
rvMC_MCR_c = zeros(row,6);
for ii = 1:row
    rTC_VNC_c_ii = rvTC_VNC_c(ii,1:3); vTC_VNC_c_ii = rvTC_VNC_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_t_ii = rvMT_MCR_t(ii,4:6);
    [rTC_MCR_c_ii,vTC_MCR_c_ii] = VNC2MCR(rTC_VNC_c_ii',vTC_VNC_c_ii',rMT_MCR_t_ii',vMT_MCR_t_ii');
    rvMC_MCR_c(ii,:) = [rTC_MCR_c_ii',vTC_MCR_c_ii'];
end

end

function [rMC_MCR_c,vMC_MCR_c] = VNC2MCR(rTC_VNC_c,vTC_VNC_c,rMT_MCR_t,vMT_MCR_t)
hT = cross(rMT_MCR_t,vMT_MCR_t);
hT_norm = norm(hT);
rMT_norm = norm(rMT_MCR_t);
vMT_norm = norm(vMT_MCR_t);
omegaT = hT/rMT_norm^2;

%VNC����ϵ����
iVNC = vMT_MCR_t/vMT_norm;
jVNC = hT/hT_norm;
kVNC = cross(iVNC,jVNC);

%����ת��
TM_VNC2MCR = [iVNC jVNC kVNC];
rTC_MCR = TM_VNC2MCR*rTC_VNC_c;
rMC_MCR_c = rTC_MCR;
% rMC_MCR_c = rMT_MCR_t+rTC_MCR;

vTC_MCR = TM_VNC2MCR*vTC_VNC_c+cross(omegaT,rTC_MCR);
vMC_MCR_c = vTC_MCR;
% vMC_MCR_c = vMT_MCR_t+vTC_MCR;
end

