function show_message(aux,locs)
fprintf('���μ���Ĺ����Ϊ: %s\n',aux.orb_type)
fprintf('���㷽��Ϊ�� �ಽ��� + α�������� + �ֲ����\n')

if strcmp(aux.orb_type , 'Lyapunov') 
    fprintf('ƽ����Ϊ: %s.\n',aux.Lyapunov_type)
end
if strcmp(aux.orb_type , 'Halo') 
    fprintf('ƽ����Ϊ: %s.\n',aux.Halo_type)
     fprintf('�������Ϊ: %s.\n',aux.Halo_class)
end
if strcmp(aux.orb_type , 'TriangleOrb') 
    fprintf('ƽ����Ϊ: %s.\n',aux.TriangleL_type)
     fprintf('�������Ϊ: %s.\n',aux.TriangleL_class)
end
if strcmp(aux.orb_type , 'ResoOrb') 
 
     fprintf('�������Ϊ: %s.\n',aux.resoOrb_type)
end
fprintf('ԭʼ�������No.%0.0f�ֲ棡\n' , aux.burfLoc)
if  nargin == 2
    if length(locs)>=2
        
        for ii = 2:length(locs)
            fprintf('3D���������µķֲ�㣬No.%0.0f��\n' , locs(ii))
        end
    else
        fprintf('3D�����û�г����µķֲ��. �ֲ����Ϊ��No.%0.0f��\n' , locs)
    end
   
end
 fprintf('--------------------�������---------------------\n')

