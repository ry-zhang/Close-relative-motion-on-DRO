function show_message(aux,locs)
fprintf('本次计算的轨道族为: %s\n',aux.orb_type)
fprintf('计算方法为： 多步打靶 + 伪弧长延拓 + 分岔分析\n')

if strcmp(aux.orb_type , 'Lyapunov') 
    fprintf('平动点为: %s.\n',aux.Lyapunov_type)
end
if strcmp(aux.orb_type , 'Halo') 
    fprintf('平动点为: %s.\n',aux.Halo_type)
     fprintf('轨道类型为: %s.\n',aux.Halo_class)
end
if strcmp(aux.orb_type , 'TriangleOrb') 
    fprintf('平动点为: %s.\n',aux.TriangleL_type)
     fprintf('轨道类型为: %s.\n',aux.TriangleL_class)
end
if strcmp(aux.orb_type , 'ResoOrb') 
 
     fprintf('轨道类型为: %s.\n',aux.resoOrb_type)
end
fprintf('原始轨道族在No.%0.0f分叉！\n' , aux.burfLoc)
if  nargin == 2
    if length(locs)>=2
        
        for ii = 2:length(locs)
            fprintf('3D轨道族出现新的分叉点，No.%0.0f！\n' , locs(ii))
        end
    else
        fprintf('3D轨道族没有出现新的分叉点. 分叉点仍为，No.%0.0f！\n' , locs)
    end
   
end
 fprintf('--------------------计算完毕---------------------\n')

