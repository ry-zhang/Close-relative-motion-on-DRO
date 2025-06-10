function [rvTC_VVLH_c] = T_MCR2VVLH(rvTC_MCR_c,rvMT_MCR_t)
%�����ٺ����������λ��rTC�����ģ�����ԭ�㣩��תϵMCR�еķ�������ת����VVLHϵ��
% v1 2020/8/20 fhl
% v2 2021/5/28 YCH
%���룺
%    ���ٺ��������λ����MCR�еķ�������rTC_MCR
%    Ŀ�꺽����λ���ٶ�rMT_MCR,vMT_MCR
%�����
%    ���ٺ��������λ����VVLH�еķ�������rTC_VVLH
%�������ϵVVLH
%    -ԭ��λ��Ŀ�꺽��������
%    -y������Ŀ�꺽�����������������Ľ��ٶȵĸ�����
%    -z������Ŀ�꺽����ָ����������
%������תϵMCR
%    -ԭ��λ������
%    -x��������ָ�������
%    -z��Ϊ˲ʱ��������ڵ���ĽǶ�������

row = size(rvTC_MCR_c,1);
rvTC_VVLH_c = zeros(row,6);
for ii = 1:row
    rTC_MCR_c_ii = rvTC_MCR_c(ii,1:3); vTC_MCR_ii = rvTC_MCR_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_ii = rvMT_MCR_t(ii,4:6);
    [rTC_VVLH_c_ii,vTC_VVLH_c_ii] = MCR2VVLH(rTC_MCR_c_ii',vTC_MCR_ii',rMT_MCR_t_ii',vMT_MCR_ii');
    rvTC_VVLH_c(ii,:) = [rTC_VVLH_c_ii',vTC_VVLH_c_ii'];
end

end

function [rTC_VVLH_c,vTC_VVLH_c] = MCR2VVLH(rTC_MCR_c,vTC_MCR_c,rMT_MCR_t,vMT_MCR_t)
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
TM_MCR2VVLH = TM_VVLH2MCR';

rTC_VVLH_c = TM_MCR2VVLH*rTC_MCR_c;

vTC_VVLH_c = TM_MCR2VVLH*(vTC_MCR_c-cross(omegaT,rTC_MCR_c));
end

