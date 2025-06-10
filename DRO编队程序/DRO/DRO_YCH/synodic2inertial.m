function rv_I_all = synodic2inertial(rv_S_all,t_all)
%% 会合坐标系至惯性系的转换
length_t = length(t_all);
rv_I_all = zeros(size(rv_S_all));
for ii_loop = 1:length_t
    theta = t_all(ii_loop);
    % 考虑平面内的转换
    if any(size(rv_S_all,1) == [2,4])
        M = [cos(theta),-sin(theta);
            sin(theta),cos(theta)];
        Mdot = [-sin(theta),-cos(theta);
                cos(theta),-sin(theta)];
        r_S = rv_S_all(1:2,ii_loop); 
        r_I = M*r_S;
        if size(rv_S_all,1)==2
            rv_I_all(:,ii_loop) = r_I;
        else
            v_S = rv_S_all(3:4,ii_loop);
            v_I = M*v_S + Mdot*r_S;
            rv_I_all(:,ii_loop) = [r_I; v_I];
        end
    else % 考虑平面外的转换
        M = [cos(theta),-sin(theta),0;
            sin(theta),cos(theta),0;
            0,0,1];
        Mdot = [-sin(theta),-cos(theta),0;
                cos(theta),-sin(theta),0;
                0,0,0];
        r_S = rv_S_all(1:3,ii_loop); 
        r_I = M*r_S;
        if size(rv_S_all,1)==3
            rv_I_all(:,ii_loop) = r_I;
        else
            v_S = rv_S_all(4:6,ii_loop);
            v_I = M*v_S + Mdot*r_S;
            rv_I_all(:,ii_loop) = [r_I; v_I];
        end
    end
end

