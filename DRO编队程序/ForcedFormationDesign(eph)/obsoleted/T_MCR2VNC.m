function [rvTC_VNC_c] = T_MCR2VNC(rvTC_MCR_c,rvMT_MCR_t)
%�����ٺ����������λ��rTC�����ģ�����ԭ�㣩��תϵMCR�еķ�������ת����VNCϵ��
% v1 2020/8/20 fhl
% v2 2021/6/1 YCH

%���룺
%    ���ٺ��������λ����MCR�еķ�������rTC_MCR
%    Ŀ�꺽����λ���ٶ�rMT_MCR,vMT_MCR
%�����
%    ���ٺ��������λ����VNC�еķ�������rTC_VNC

%�������ϵVNC
%    -ԭ��λ��Ŀ�꺽��������
%    -x��������������ָ��Ŀ�꺽����
%    -z������Ŀ�꺽�����������������Ľ��ٶȷ���
%������תϵMCR
%    -ԭ��λ������
%    -x��������ָ�������
%    -z��Ϊ˲ʱ��������ڵ���ĽǶ�������

row = size(rvTC_MCR_c,1);
rvTC_VNC_c = zeros(row,6);
for ii = 1:row
    rTC_MCR_c_ii = rvTC_MCR_c(ii,1:3); vTC_MCR_ii = rvTC_MCR_c(ii,4:6);
    rMT_MCR_t_ii = rvMT_MCR_t(ii,1:3); vMT_MCR_ii = rvMT_MCR_t(ii,4:6);
    [rTC_VNC_c_ii,vTC_VNC_c_ii] = MCR2VNC(rTC_MCR_c_ii',vTC_MCR_ii',rMT_MCR_t_ii',vMT_MCR_ii');
    rvTC_VNC_c(ii,:) = [rTC_VNC_c_ii',vTC_VNC_c_ii'];
end

end

function [rTC_VNC_c,vTC_VNC_c] = MCR2VNC(rTC_MCR_c,vTC_MCR_c,rMT_MCR_t,vMT_MCR_t)
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
TM_MCR2VNC = TM_VNC2MCR';

rTC_VNC_c = TM_MCR2VNC*rTC_MCR_c;

vTC_VNC_c = TM_MCR2VNC*(vTC_MCR_c-cross(omegaT,rTC_MCR_c));
end

