function rv_S_all = inertial2synodic(rv_I_all,t_all)
%% 会合坐标系至惯性系的转换
%% 会合坐标系至惯性系的转换
length_t = length(t_all);
rv_S_all = zeros(size(rv_I_all));
for ii_loop = 1:length_t
    theta = t_all(ii_loop);
    % 考虑平面内的转换
    if any(size(rv_I_all,1) == [2,4])
        M = [cos(theta),-sin(theta);
            sin(theta),cos(theta)];
        Mdot = [-sin(theta),-cos(theta);
                cos(theta),-sin(theta)];
        r_I = rv_I_all(1:2,ii_loop); 
        r_S = M'*r_I;
        if size(rv_I_all,1)==2
            rv_S_all(:,ii_loop) = r_S;
        else
            v_I = rv_I_all(3:4,ii_loop);
            v_S = M'*v_I + Mdot'*r_I;
            rv_S_all(:,ii_loop) = [r_S; v_S];
        end
    else % 考虑平面外的转换
        M = [cos(theta),-sin(theta),0;
            sin(theta),cos(theta),0;
            0,0,1];
        Mdot = [-sin(theta),-cos(theta),0;
                cos(theta),-sin(theta),0;
                0,0,0];
        r_I = rv_I_all(1:3,ii_loop); 
        r_S = M'*r_I;
        if size(rv_I_all,1)==3
            rv_S_all(:,ii_loop) = r_S;
        else
            v_I = rv_I_all(4:6,ii_loop);
            v_S = M'*v_I + Mdot'*r_I;
            rv_S_all(:,ii_loop) = [r_S; v_S];
        end
    end
end
