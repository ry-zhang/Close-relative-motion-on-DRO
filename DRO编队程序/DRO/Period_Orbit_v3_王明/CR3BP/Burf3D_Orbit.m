function dataSave = Burf3D_Orbit(nLoop,currBurf , aux)
% dataSave = zeros(nLoop, 9);
im_Loop = 0;
for iLoop = 1 : nLoop
    % 预测
    xjp1 = aux.xj + aux.delta_xj * aux.delta_s;
    % 修正
    err = inf;
    count = 1;
    iterMax = 20;
    while err > 1e-10 && count < iterMax
        % 约束 / 约束梯度
        [~, cc_eq, ~, Gc_eq_] = constr(xjp1, aux);
        Gc_eq = Gc_eq_';
        
        % 计算误差
        h = cc_eq;
        
        err = max(abs(h));
        Jac = [Gc_eq];
        
        dx =  -inv(Jac'*Jac)*Jac'*h;
        xjp1 = xjp1 +dx;
        
        count = count + 1;
    end
    if count < iterMax
        fprintf('No: %0.0f 收敛！\n' , iLoop)
    else
        fprintf('No: %0.0f 不收敛！\n' , iLoop)
    end
    
     % xj 和 delta_xj更新
    aux.xj = xjp1;
    vOld = aux.delta_xj ;
   
    
    [~,S,V] = svd(Gc_eq(1:end-1, :));
    
    svd2ndLast = S(end, end-1);                       % 最小奇异值
    
    aux.delta_xj = V(:, end);                                % 零空间向量
      % burf searching
%     if svd2ndLast < aux.svdBurfTol 
%         vTemp =   V(:, end-1);
%         if aux.UseSymmetric_IO == 1
%             if  abs(dot(vTemp(1:3), [0,1,0]'))/norm(vTemp(1:3))>0.9
%                 fprintf('No: %0.0f 潜在3d_z分叉点！\n' , iLoop)
%             end
%         else
%             if  abs(dot(vTemp(1:6),[0,0,0, 0,0,1]'))/norm(vTemp(1:6))>0.9 || abs(dot(vTemp(1:6), [0,0,1, 0,0,0]'))/norm(vTemp(1:6))>0.9
%                 fprintf('No: %0.0f 潜在3d_zd分叉点！\n' , iLoop)
%             end
%         end
%     end
%     if iLoop==197
%          currBurf = 0;
%     end
    if svd2ndLast<aux.svdBurfTol  && aux.computeBurf &&  ~currBurf  && iLoop>=aux.burfLoc
        aux.delta_xj = V(:, end-1);
        aux.delta_s = aux.burf_delta_s;
        currBurf = 1;
    end
    
    if dot(vOld, aux.delta_xj ) <0
        aux.delta_s  = -1*aux.delta_s ;
    end
     %【测试】
     figure(2);
%      set(gcf,'Resize','off');
     % [~ , ~] = crtbpFlow3D(xjp1(1:6) , xjp1(end) , 'b' , 1 , aux);
     if aux.UseSymmetric_IO == 1
         x0 = [xjp1(1),0,xjp1(2),0,xjp1(3),0]';
         TT = 2*xjp1(end);
     else
         x0 = xjp1(1:6);
         TT = xjp1(end);
     end
     options = odeset('Reltol' , aux.tol , 'AbsTol' , aux.tol,'event',@event_impact);
     [tt_rot , xx_rot,te,xe,ie] = ode113(@crtbpEqm3D , [0 , TT] , x0 , options , aux);
     if isempty(ie)
         plot3(xx_rot(: , 1) , xx_rot(: , 2) , xx_rot(: , 3) , 'b' , 'linewidth' , 0.2);
         if iLoop > aux.burfLoc 
             plot3(xx_rot(: , 1) , xx_rot(: , 2) , xx_rot(: , 3) , 'r' , 'linewidth' , 0.2);
             
         end
   
         
         % 调整坐标系范围
         xlu_rot(iLoop,:) = [min(xx_rot(:,1:3)),max(xx_rot(:,1:3))];
         if iLoop == nLoop
             set_axis(xlu_rot);
         end
     else
         im_Loop = iLoop;
         set_axis(xlu_rot);
         fprintf('【注意！！】延拓轨道族撞击天体!!! \n');
        
     end
     % ------------------------- plot in inertial frame  -------------------------
    
     figure(1);
     [xx_ini , tt_ini] = crtbpSynodic2P1centered3D(xx_rot , tt_rot , aux.mu);
     plot3(xx_ini(: , 1) , xx_ini(: , 2) , xx_ini(: , 3) , 'b' , 'linewidth' , 0.2);
      if iLoop > aux.burfLoc
         plot3(xx_ini(: , 1) , xx_ini(: , 2) , xx_ini(: , 3) , 'r' , 'linewidth' , 0.2);
      end
      xlu_ini(iLoop,:) = [min(xx_ini(:,1:3)),max(xx_ini(:,1:3))];
      dataSave(iLoop,:) = [x0', TT, crtbpJacobi3D(x0' , aux.mu), svd2ndLast];
      % 调整坐标系范围
      if iLoop == nLoop ||  iLoop == im_Loop
          set_axis(xlu_ini);
          break;
      end
     
end
%写文件 
% filePath = [pwd , '\Orbit_data\' , '3DResoOrb1_2'];
% fid = fopen(filePath , 'a');
% fprintf(fid , '%s   %s   %s   %s   %s   %s   %s  %s\n', 'X(LU)','Y(LU)','Z(LU)','Vx(VU)','Vy(VU)','Vz(VU)','T(TU)','C');
% for ii = aux.burfLoc:size(dataSave,1)
%     fprintf(fid , '%5.15f %5.15f %5.15f %5.15f %5.15f %5.15f %5.15f %5.15f \n' ,  dataSave(ii,1:8)');
%     
% end
% fclose('all');







[~,locs] = findpeaks(-dataSave(:,9));
locs = locs(locs.*(dataSave(locs,9)<aux.svdBurfTol)>0);
show_message(aux,locs);
 



%%
figure(3);hold on
plot(dataSave(:,7) , 'r');
grid on
ylabel('P')
xlabel('# of NumCont')

figure(4);hold on
plot(dataSave(:,8) , 'r');
grid on
ylabel('C')
xlabel('# of NumCont')

figure(5);hold on
plot(dataSave(:,9) , 'r'); hold on
plot(locs, dataSave(locs,9) , 'go'); hold on
ylabel('svd')
grid on
xlabel('# of NumCont')
end

