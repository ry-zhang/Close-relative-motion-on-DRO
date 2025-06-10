function [locs,dataSave] = Origin_Orbit_family(nLoop, aux)
%
 im_Loop = 0;
% dataSave = zeros(nLoop, 9);
for iLoop = 1 : nLoop
    % Ԥ��
    xjp1 = aux.xj + aux.delta_xj * aux.delta_s;
    % ����
    err = inf;
    count = 1;
    iterMax = 20;
    while err > 1e-10 && count < iterMax
        % Լ�� / Լ���ݶ�
        [~, cc_eq, ~, Gc_eq_] = constr(xjp1, aux);
        Gc_eq = Gc_eq_';
        
        % �������
        h = cc_eq;
        
        err = max(abs(h));
        Jac = [Gc_eq];
        
        dx =  -inv(Jac'*Jac)*Jac'*h;
        xjp1 = xjp1 +dx;
        
        count = count + 1;
    end
    if count < iterMax
        fprintf('No: %0.0f ������\n' , iLoop)
    else
        fprintf('No: %0.0f ��������\n' , iLoop)
    end
    
     % xj �� delta_xj����
    aux.xj = xjp1;
    vOld = aux.delta_xj ;
   
    
    [~,S,V] = svd(Gc_eq(1:end-1, :));
    
    svd2ndLast = S(end, end-1);                       % ��С����ֵ
    
    aux.delta_xj = V(:, end);                                % ��ռ�����
      % burf searching
    if svd2ndLast < aux.svdBurfTol 
        vTemp =   V(:, end-1);
        if aux.UseSymmetric_IO == 1
            if  abs(dot(vTemp(1:3), [0,1,0]'))/norm(vTemp(1:3))>0.9
                fprintf('No: %0.0f Ǳ��3d_z�ֲ�㣡\n' , iLoop)
            end
        else
            if  abs(dot(vTemp(1:6),[0,0,0, 0,0,1]'))/norm(vTemp(1:6))>0.9 || abs(dot(vTemp(1:6), [0,0,1, 0,0,0]'))/norm(vTemp(1:6))>0.9
                fprintf('No: %0.0f Ǳ��3d_zd�ֲ�㣡\n' , iLoop)
            end
        end
    end
    
    if dot(vOld, aux.delta_xj ) <0
        aux.delta_s  = -1*aux.delta_s ;
    end
     %�����ԡ�
     figure(2);
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
         xlu_rot(iLoop,:) = [min(xx_rot(:,1:3)),max(xx_rot(:,1:3))];
   
         % ��������ϵ��Χ
         
     else
         im_Loop = iLoop;
         fprintf('��ע�⣡�������ع����ײ������!!! \n');
         set_axis(xlu_rot);
        
     end
     if iLoop == nLoop 
         set_axis(xlu_rot);
     end
      
     % ------------------------- plot in inertial frame  -------------------------
     figure(1);
     [xx_ini , tt_ini] = crtbpSynodic2P1centered3D(xx_rot , tt_rot , aux.mu);
     plot3(xx_ini(: , 1), xx_ini(: , 2) , xx_ini(: , 3) , 'b' , 'linewidth' , 0.2);
     xlu_ini(iLoop,:) = [min(xx_ini(:,1:3)),max(xx_ini(:,1:3))];
     dataSave(iLoop,:) = [x0', TT, crtbpJacobi3D(x0' , aux.mu), svd2ndLast];
     % ��������ϵ��Χ
     if iLoop == nLoop || iLoop == im_Loop
         set_axis(xlu_ini);
         break;
     end
    
end
%д�ļ� 
filePath = [pwd , '\Orbit_data\' , 'RO_3_4_2D.txt'];
fid = fopen(filePath , 'a');
fprintf(fid , '%s   %s   %s   %s   %s   %s   %s  %s\n', 'X(LU)','Y(LU)','Z(LU)','Vx(VU)','Vy(VU)','Vz(VU)','T(TU)','C');
for ii = 1:size(dataSave,1)
    fprintf(fid , '%5.15f %5.15f %5.15f %5.15f %5.15f %5.15f %5.15f %5.15f \n' ,  dataSave(ii,1:8)');
    
end
fclose('all');
[~,locs] = findpeaks(-dataSave(:,9));
locs = locs(locs.*(dataSave(locs,9)<aux.svdBurfTol)>0);
fprintf('No: %0.0f Ǳ�ڷֲ�㣡\n' , locs)
aux.burfLoc = locs;
show_message(aux);

%%
figure(3);
plot(dataSave(:,7) , 'b');
grid on
ylabel('P')
xlabel('# of NumCont')

figure(4);
plot(dataSave(:,8) , 'b');
grid on
ylabel('C')
xlabel('# of NumCont')

figure(5);
plot(dataSave(:,9) , 'b'); hold on
plot(locs, dataSave(locs,9) , 'ro'); hold on
ylabel('svd')
grid on
xlabel('# of NumCont')

end

