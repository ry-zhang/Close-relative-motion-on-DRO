
function [yy , aux] = getIC_periOrb(aux)
%
% 单步打靶变多步打靶
% 
% 例子：
% n节点周期轨道：
% A---------B
% |             |
% |             |
% C---------D
% 变量数：6 * n(rv) + 1(T) = 25
% 
% 作者：张晨
% 2021年5月30日
%%%%%%%%%%%%%%%%%%%%%%%%

% 周期轨道节点数
n = aux.node_n;
mu = aux.mu;
% aux.Halopseudo = 0;                                         % Halo对伪弧步数敏感，所以单独设置
% 构造积分初值

% 建议打靶方式

if strcmp(aux.orb_type , 'DRO') || strcmp(aux.orb_type , 'Lyapunov') || strcmp(aux.orb_type , 'DPO') || strcmp(aux.orb_type , 'LoPO')
    if aux.UseSymmetric_IO == 0
        fprintf('--------------------------------------------------- \n');
        prompt = ('建议您使用轨道对称性进行校正，以提高效率，请输入字符串: ''yes''或数字：1更改 \n');
        
        str = input(prompt);
        
        
        if strcmp(str , 'yes') || strcmp(str , 'YES') || strcmp(str , 'Yes') || str==1
            aux.UseSymmetric_IO = 1;
        else
            fprintf('--------------------------------------------------- \n');
            fprintf('本次已自动帮您更改aux.UseSymmetric_IO = 1! \n');
            pause(2);
        end
        
    end
elseif strcmp(aux.orb_type , 'ResoOrb') || strcmp(aux.orb_type , 'TriangleOrb')
    if aux.UseSymmetric_IO == 1
        fprintf('--------------------------------------------------- \n');
        prompt = ('由于轨道非对称，aux.UseSymmetric_IO 应为 0 ，请输入字符串: ''yes''或数字：1更改 \n');
        str = input(prompt);
        if strcmp(str , 'yes') || strcmp(str , 'YES') || strcmp(str , 'Yes') || str==1
            aux.UseSymmetric_IO = 0;
        else
             fprintf('--------------------------------------------------- \n');
             fprintf('已自动帮您更改aux.UseSymmetric_IO = 0! \n');
             pause(2);
             aux.UseSymmetric_IO = 0;
        end
    end
elseif strcmp(aux.orb_type , 'Halo')
    if aux.UseSymmetric_IO == 0 
         fprintf('--------------------------------------------------- \n');
         fprintf('检查程序确定aux.UseSymmetric_IO = 1，本次已自动帮您更改! \n');
         pause(1);
         aux.UseSymmetric_IO = 1;
    end
         if  aux.Halopseudo == 0
             fprintf('--------------------------------------------------- \n');
             fprintf('检查程序确定aux.Halopseudo == 1，本次已自动帮您更改! \n');
             pause(2);
             aux.Halopseudo = 1;
        end
end
  
if strcmp(aux.orb_type , 'DRO')
  
    x0 = aux.DRO_xx;                                    % 轨道初值
    
    T = aux.DRO_P;
    
elseif strcmp(aux.orb_type , 'Halo')
    %     aux.Halopseudo = 1;
%     aux.Halopseudo = 1;
    % Richardson三阶近似解析解
    Az = 0.01;                                                   % 法向振幅
    [L1, L2, ~, ~, ~] = Li (aux.mu);                      % 平动点坐标
    if strcmp(aux.Halo_type , 'L1')
        x0 = FirstGuessFixedAz3D (Az, L1, mu, aux.Halo_type, aux.Halo_class);
    elseif strcmp(aux.Halo_type , 'L2')
        x0 = FirstGuessFixedAz3D (Az, L2, mu, aux.Halo_type, aux.Halo_class);
    end
     % 利用过y轴检测估计估计周期
     options = odeset('Reltol', 1e-12 , 'AbsTol' , 1e-12,'event',@ev_orbit3D);
     [tt , xx,te,xe,ie] = ode113(@crtbpEqm3D, [0 , pi] , x0 , options , aux);
     if isempty(ie)
         warning('Unable to compute first guess solution') ;
         return;
     end
     
     T = 2*te;
   
elseif strcmp(aux.orb_type , 'Lyapunov')
     if strcmp(aux.Lyapunov_type , 'L1')
         x0 = aux.Lyapunov_xx(1,:);
         T = aux.Lyapunov_P(1);
     elseif strcmp(aux.Lyapunov_type , 'L2')
           x0 = aux.Lyapunov_xx(2,:);
           T = aux.Lyapunov_P(2);
     elseif strcmp(aux.Lyapunov_type , 'L3')
         x0 = aux.Lyapunov_xx(3,:);
         T = aux.Lyapunov_P(3);
     end
elseif strcmp(aux.orb_type , 'ResoOrb')
%     aux.node_n = 20;
    
     [yy , aux] = getIC_resoOrb(aux);
     return;
elseif strcmp(aux.orb_type , 'TriangleOrb')
    if strcmp(aux.TriangleL_type , 'L4')
        if strcmp(aux.TriangleL_class , 'Short_P')
            x0 = aux.TriangleL_xx(1,:);
            T = aux.TriangleL_P(1);
        elseif strcmp(aux.TriangleL_class , 'Long_P')
            x0 = aux.TriangleL_xx(2,:);
            T = aux.TriangleL_P(2);
        end
    elseif strcmp(aux.TriangleL_type , 'L5')
         if strcmp(aux.TriangleL_class , 'Short_P')
             x0 = aux.TriangleL_xx(3,:);
             T = aux.TriangleL_P(3);
        elseif strcmp(aux.TriangleL_class , 'Long_P')
            x0 = aux.TriangleL_xx(4,:);
            T = aux.TriangleL_P(4);
        end
    end
elseif strcmp(aux.orb_type , 'DPO')
    x0 = aux.DPO_xx;                                    % 轨道初值
    
    T = aux.DPO_P;
elseif strcmp(aux.orb_type , 'LoPO')
    x0 = aux.LoPO_xx;                                    % 轨道初值
    
    T = aux.LoPO_P; 
else
    fprintf('wrong!');
end
if aux.UseSymmetric_IO == 1
    aux.periOrb_P = T/2;
else
    aux.periOrb_P = T;
end
% 数值积分
options = odeset('Reltol', 1e-12 , 'AbsTol' , 1e-12);
[tt , xx] = ode113(@crtbpEqm3D, [0 , aux.periOrb_P] , x0 , options , aux);

if aux.Halopseudo == 1
    n = 1;
end
% 三次样条插值
tt_interp = linspace(0 , aux.periOrb_P , n + 1)';
xx_interp = interp1(tt , xx , tt_interp' , 'spline') ;

% 构造优化变量
yy = zeros(6 * n + 1 ,1) ;
for iLoop = 1 : n
   yy(6 * (iLoop - 1) + 1 : 6 * iLoop) = xx_interp(iLoop , :)' ;
end

% 周期轨道周期5
yy(6 * n + 1) = aux.periOrb_P;

% 固定共振轨道第一个点的y值！
aux.yFix = yy(2);

if aux.UseSymmetric_IO == 1
    yy([2,4,6]) = [];                                               % 如果利用对称性，自由变量去掉y, x_dot, z_dot
end
% ----------------------- 【测试】-----------------------
% figure(1); hold on; grid on;
% crtbpMarkEM;
% aux.h_g = plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'g' , 'linewidth' , 1);
% plot3(xx_interp(1:end - 1 , 1) , xx_interp(1:end - 1 , 2) , xx_interp(1:end - 1 , 3) , '.r' , 'linewidth' ,2,'Markersize',15);

end
