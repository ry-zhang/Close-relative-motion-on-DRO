function [rvTC_LVLH_c] = T_MCR2LVLH(rvTC_MCR_c,rvMT_MCR_t)
%�����ٺ����������λ��rTC�����ģ�����ԭ�㣩��תϵMCR�еķ�������ת����LVLHϵ��
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH
% LVLH����ϵ��RIC

%���룺
%    ���ٺ��������λ����MCR�еķ�������rTC_MCR
%    Ŀ�꺽����λ���ٶ�rMT_MCR,vMT_MCR
%�����
%    ���ٺ��������λ����LVLH�еķ�������rTC_LVLH

%�������ϵLVLH
%    -ԭ��λ��Ŀ�꺽��������
%    -x��������������ָ��Ŀ�꺽����
%    -z������Ŀ�꺽�����������������Ľ��ٶȷ���
%������תϵMCR
%    -ԭ��λ������
%    -x��������ָ�������
%    -z��Ϊ˲ʱ��������ڵ���ĽǶ�������

row = size(rvTC_MCR_c,1);
rvTC_LVLH_c = zeros(row,6);
for ii = 1:row
    rTC_MCR_c_ii = rvTC_MCR_c(ii,1:3); vTC_MCR_ii = rvTC_MCR_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_ii = rvMT_MCR_t(ii,4:6);
    [rTC_LVLH_c_ii,vTC_LVLH_c_ii] = MCR2LVLH(rTC_MCR_c_ii',vTC_MCR_ii',rMT_MCR_t_ii',vMT_MCR_ii');
    rvTC_LVLH_c(ii,:) = [rTC_LVLH_c_ii',vTC_LVLH_c_ii'];
end

end

function [rTC_LVLH_c,vTC_LVLH_c] = MCR2LVLH(rTC_MCR_c,vTC_MCR_c,rMT_MCR_t,vMT_MCR_t)
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
TM_MCR2LVLH = TM_LVLH2MCR';

rTC_LVLH_c = TM_MCR2LVLH*rTC_MCR_c;

vTC_LVLH_c = TM_MCR2LVLH*(vTC_MCR_c-cross(omegaT,rTC_MCR_c));
end

