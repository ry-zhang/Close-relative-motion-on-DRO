% �Ѿ���
clear
clc
%�����ʼ��Ԫ�����գ�UTCʱ��
load('data_DROSFF_withOrbit2.mat');
load('structState.mat')
fid=fopen('jd0UTC.txt','wt');
fprintf(fid,'%f',aux_SFF.jd0UTC);
fclose(fid);
%����������
fid_A=fopen('orbit-a.txt','wt');
% fid_B=fopen('orbit-b.txt','wt');
fid_AB=fopen('deltaX-ab.txt','wt');
fid_B=fopen('orbit-b.txt','wt');
%����廬�ж�
% for i=1:3707
%     fprintf(fid_A,'%f ',aux_SFF.Coasting(1).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.Coasting(1).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.Coasting(1).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.Coasting(1).xx_j2k_AB(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.Coasting(1).xx_j2k_AB(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.Coasting(1).xx_j2k_AB(i,j)*1000-aux_SFF.Coasting(1).xx_j2k_AB(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
% %�����
% for i=1:145
%     fprintf(fid_A,'%f ',aux_SFF.transferSFF(1).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.transferSFF(1).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.transferSFF(1).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.transferSFF(1).xx_j2k_A(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.transferSFF(1).xx_j2k_B(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.transferSFF(1).xx_j2k_B(i,j)*1000-aux_SFF.transferSFF(1).xx_j2k_A(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
% %�����ĩ�� �� �̾����ɶε�ת��
% for i=1:1441
%     fprintf(fid_A,'%f ',aux_SFF.transferSFF(2).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.transferSFF(2).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.transferSFF(2).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.transferSFF(2).xx_j2k_A(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.transferSFF(2).xx_j2k_B(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.transferSFF(2).xx_j2k_B(i,j)*1000-aux_SFF.transferSFF(2).xx_j2k_A(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
%�̾����ɶ�
% for i=1:7525
%     fprintf(fid_A,'%f ',aux_SFF.naturalSFF(1).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.naturalSFF(1).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.naturalSFF(1).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.naturalSFF(1).xx_j2k_A(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.naturalSFF(1).xx_j2k_B(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.naturalSFF(1).xx_j2k_B(i,j)*1000-aux_SFF.naturalSFF(1).xx_j2k_A(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
for i=1:21601
    fprintf(fid_A,'%f ',aux_SFF.jd0UTC+structState.t_sample(i));
    fprintf(fid_B,'%f ',aux_SFF.jd0UTC+structState.t_sample(i));
    fprintf(fid_AB,'%f ',aux_SFF.jd0UTC+structState.t_sample(i));
    for j=1:6
        fprintf(fid_A,'%f ',structState.xx_MCR_target(i,j)*1000);
        fprintf(fid_B,'%f ',structState.xx_MCR_chaser(i,j)*1000);
        fprintf(fid_AB,'%f ',structState.xx_MCRLVLH_rel(i,j)*1000);
    end
    fprintf(fid_A,'\n');
    fprintf(fid_B,'\n');
    fprintf(fid_AB,'\n');
end
% %�̾����ɶ� �� �еȾ����ɶε�ת��
% for i=1:1441
%     fprintf(fid_A,'%f ',aux_SFF.transferSFF(3).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.transferSFF(3).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.transferSFF(3).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.transferSFF(3).xx_j2k_A(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.transferSFF(3).xx_j2k_B(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.transferSFF(3).xx_j2k_B(i,j)*1000-aux_SFF.transferSFF(3).xx_j2k_A(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
% %�еȾ����ɶ�
% for i=1:6949
%     fprintf(fid_A,'%f ',aux_SFF.naturalSFF(2).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.naturalSFF(2).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.naturalSFF(2).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.naturalSFF(2).xx_j2k_A(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.naturalSFF(2).xx_j2k_B(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.naturalSFF(2).xx_j2k_B(i,j)*1000-aux_SFF.naturalSFF(2).xx_j2k_A(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
% %�еȾ����ɶ� �� Զ�����ɶε�ת��
% for i=1:4321
%     fprintf(fid_A,'%f ',aux_SFF.transferSFF(4).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.transferSFF(4).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.transferSFF(4).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.transferSFF(4).xx_j2k_A(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.transferSFF(4).xx_j2k_B(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.transferSFF(4).xx_j2k_B(i,j)*1000-aux_SFF.transferSFF(4).xx_j2k_A(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
% %Զ�����ɶ�
% for i=1:7669
%     fprintf(fid_A,'%f ',aux_SFF.naturalSFF(3).t_sample(i));
%     fprintf(fid_B,'%f ',aux_SFF.naturalSFF(3).t_sample(i));
%     fprintf(fid_AB,'%f ',aux_SFF.naturalSFF(3).t_sample(i));
%     for j=1:6
%         fprintf(fid_A,'%f ',aux_SFF.naturalSFF(3).xx_j2k_A(i,j)*1000);
%         fprintf(fid_B,'%f ',aux_SFF.naturalSFF(3).xx_j2k_B(i,j)*1000);
%         fprintf(fid_AB,'%f ',aux_SFF.naturalSFF(3).xx_j2k_B(i,j)*1000-aux_SFF.naturalSFF(3).xx_j2k_A(i,j)*1000);
%     end
%     fprintf(fid_A,'\n');
%     fprintf(fid_B,'\n');
%     fprintf(fid_AB,'\n');
% end
% fclose(fid_A);
% fclose(fid_B);
% fclose(fid_AB);
% %�����������
% fid=fopen('impulse-b.txt','wt');
% %������������
% fprintf(fid,'%f ',aux_SFF.transferSFF(1).tm0);
% for i=1:3
%     fprintf(fid,'%f ',aux_SFF.transferSFF(1).dvm0_j2k(i));
% end
% fprintf(fid,'\n');
% %�����ĩ��
% fprintf(fid,'%f ',aux_SFF.transferSFF(2).tm0);
% for i=1:3
%     fprintf(fid,'%f ',aux_SFF.transferSFF(2).dvm0_j2k(i));
% end
% fprintf(fid,'\n');
% %����̾����ɶ�����
% fprintf(fid,'%f ',aux_SFF.transferSFF(2).tmf);
% for i=1:3
%     fprintf(fid,'%f ',aux_SFF.transferSFF(2).dvmf_j2k(i));
% end
% fprintf(fid,'\n');
% %�뿪�̾����ɶ�����
% fprintf(fid,'%f ',aux_SFF.transferSFF(3).tm0);
% for i=1:3
%     fprintf(fid,'%f ',aux_SFF.transferSFF(3).dvm0_j2k(i));
% end
% fprintf(fid,'\n');
% %�����еȾ����ɶ�����
% fprintf(fid,'%f ',aux_SFF.transferSFF(3).tmf);
% for i=1:3
%     fprintf(fid,'%f ',aux_SFF.transferSFF(3).dvmf_j2k(i));
% end
% fprintf(fid,'\n');
% %�뿪�еȾ����ɶ�����
% fprintf(fid,'%f ',aux_SFF.transferSFF(4).tm0);
% for i=1:3
%     fprintf(fid,'%f ',aux_SFF.transferSFF(4).dvm0_j2k(i));
% end
% fprintf(fid,'\n');
% %����Զ�����ɶ�����
% fprintf(fid,'%f ',aux_SFF.transferSFF(4).tmf);
% for i=1:3
%     fprintf(fid,'%f ',aux_SFF.transferSFF(4).dvmf_j2k(i));
% end
% fprintf(fid,'\n');
%�������
fclose(fid);

%% load���ݣ�������
load('orbit-a.txt');
plot(orbit_a(:,1),orbit_a(:,2))
xlim(aux_SFF.transferSFF(1).tm0+[-180,180])
% xlim(aux_SFF.transferSFF(2).tm0+[-180,180])
% xlim(aux_SFF.transferSFF(3).tm0+[-180,180])
% xlim(aux_SFF.transferSFF(4).tm0+[-180,180])