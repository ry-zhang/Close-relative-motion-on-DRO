function [yy,aux] = getIC_periOrb(aux)
%
% 获得周期轨道初值
%
% 作者：wm
% 2021年11月3日
%%%%%%%%%%%%%%%%%%%%%%%%
% 读数据
addpath('[20211104]_crtbpFamily_v2');
FileName = aux.periOrb_Modi +".txt";
fid = fopen(FileName , 'r');
tline = fgets(fid);
tline = fgets(fid);
count = 1;
while tline ~= -1
    xxNow = textscan(tline,'%f');
    xxNow = cell2mat(xxNow);
    X0Mtx(count,:) = xxNow';
    tline = fgets(fid);
    count = count +1 ;
end
fclose('all');

% 该族轨道周期范围
Pmin = min(X0Mtx(: , end - 1));
Pmax = max(X0Mtx(: , end - 1));
% 该族轨道能量范围
Cmin = min(X0Mtx(: , end));
Cmax = max(X0Mtx(: , end));
fprintf('轨道族类型：%s \n' , aux.periOrb_Modi)
fprintf('轨道族周期：[%0.4f , %0.4f] \n' , Pmin , Pmax )
fprintf('轨道族能量：[%0.4f , %0.4f] \n' , Cmin , Cmax)
fprintf('--------------------------------------------------- \n');
% 月球周期倍数
periOrb_m = aux.periOrb_m;

% 周期轨道一圈节点数
periOrb_n = aux.periOrb_n;
% 判断是否共振轨道
if strcmp(aux.periOrb_Modi(1:2) , 'RO')
    num_str = regexp(aux.periOrb_Modi,'\d*\.?\d*','match');
    num = str2double(num_str);
    if isempty(aux.periOrb_T)
        
        
        periOrb_p = num(1);
        periOrb_q = num(2);
        TT = 2*pi*periOrb_q;
        
    else
        TT = str2double(aux.periOrb_T);
        periOrb_p =num(1);
    end
    
    Orbit_data = X0Mtx;
    if TT>max(X0Mtx(:,7)) || TT<min(X0Mtx(:,7))
        fprintf('--------------------------------------------------- \n');
        fprintf('当前周期对应轨道不存在! \n');
        fprintf('建议您输入周期范围：[%0.4f , %0.4f],' , Pmin, Pmax )
        prompt = ('请重新输入轨道周期： \n');
        num = input(prompt);
        if num<=Pmax && num>=Pmin
            TT = num;
        else
            fprintf('输入错误，程序终止! \n');
            return;
        end
    end
    % 对插值进行打靶校正
    fprintf(' ------- 三体初值校正 --------- \n')
    x0 =  interp1(Orbit_data(:,7) , Orbit_data(:,1:6) , TT , 'spline');  %    轨道初始位置
    x0 = Refine_x0(x0, TT, aux);
    options = odeset('RelTol',1e-13,'AbsTol',1e-16);
    [tt, xx]  = ode113(@crtbpEqm3D, [0, TT*periOrb_m], x0 , options, aux);
    figure (100)
    plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'b');
    % 三次样条插值
    tt_interp = linspace(0 , TT*periOrb_m , periOrb_n*periOrb_m*periOrb_p)';
    xx_interp = interp1(tt , xx , tt_interp' , 'spline') ;
%     fprintf('初始太阳相角：%0.0f° \n' , rad2deg(aux.Sun_angle))
    fprintf('轨道转数：%0.0f，初始周期 %0.4f\n' , periOrb_m*periOrb_p,TT*periOrb_m)
else
   
    
    if isempty(aux.periOrb_T)
        % 判断 p,q是否公约
        [~, gcd] = is_coprime(aux.resoOrb_p  , aux.resoOrb_q);
        
        % 共振比，p : q
        periOrb_p = aux.periOrb_p/gcd ;
        periOrb_q  = aux.periOrb_q/gcd ;
        TT = 2*pi*periOrb_q/periOrb_p;                             % 单圈轨道周期
        fprintf('%0.0f ：%0.0f 单圈轨道周期：%0.4f \n' , periOrb_p , periOrb_q,TT)
%         TT = 1/abs(aux.oms)*2*pi*periOrb_q/periOrb_p;
%         % 判断该共振比轨道是否存在
%         fprintf('考虑太阳周期共振后，%0.0f ：%0.0f 单圈轨道周期：%0.4f \n' , periOrb_p , periOrb_q,TT)
    else
        periOrb_p = 1;
        %         periOrb_q = 1;
        TT = str2double(aux.periOrb_T);
        fprintf('%0.0f 圈轨道周期：%0.4f \n' , periOrb_p ,TT)
        %      TT = 1/abs(aux.oms)*TT;
        %     fprintf('考虑四体下太阳周期共振后，%0.0f ：%0.0f 单圈轨道周期修正为：%0.4f \n' , periOrb_p , periOrb_q,TT)
        
    end
    
    
    % 初始时刻太阳相角
    
    
    % 匹配单圈周期
    
    % 注意，这里与朔望日共振
    
    
    if TT>max(X0Mtx(:,7)) || TT<min(X0Mtx(:,7))
        periOrb_p = 1;                                              % 不涉及共振比
%         periOrb_q = 1;
        fprintf('--------------------------------------------------- \n');
        fprintf('当前周期对应轨道不存在! \n');
        fprintf('建议您输入周期范围：[%0.4f , %0.4f],' , Pmin, Pmax )
        prompt = ('请重新输入轨道周期： \n');
        num = input(prompt);
        if num<=Pmax && num>=Pmin
            TT = num;
        else
            fprintf('输入错误，程序终止! \n');
            return;
        end
    end
%     Sun_angle0 = aux.Sun_angle;
%     fprintf('初始太阳相角：%0.0f° \n' , rad2deg(Sun_angle0))
    fprintf('轨道转数：%0.0f，初始周期 %0.4f\n' , periOrb_m*periOrb_p,TT*periOrb_p*periOrb_m)
    Orbit_data = X0Mtx;
    % % 判断待修正轨道类型
    % if strcmp(aux.orb_type , 'DRO')
    %     Orbit_data = aux.DRO;
    % end
    
    x0 =  interp1(Orbit_data(:,7) , Orbit_data(:,1:6) , TT , 'spline');  %    轨道初始位置
  
    % 对插值进行打靶校正
    fprintf(' ------- 三体初值校正 --------- \n')
    x0 = Refine_x0(x0, TT, aux);
    Tspan = [0, TT*periOrb_p*periOrb_m];
    options = odeset('RelTol',1e-13,'AbsTol',1e-16);
    [tt, xx]  = ode113(@crtbpEqm3D, Tspan, x0 , options, aux);
    
    % 三次样条插值
    tt_interp = linspace(0 , TT*periOrb_p*periOrb_m , periOrb_n*periOrb_m*periOrb_p)';
    xx_interp = interp1(tt , xx , tt_interp' , 'spline') ;
    
    
    
    % 构造优化变量
end
% 【测试】
figure(1); hold on; grid on;
crtbpMarkEM;
plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'b');
plot3(xx_interp(: , 1) , xx_interp(: , 2) , xx_interp(: , 3) , 'bo');
    
% xx_interp = xx_interp';
n = periOrb_n*periOrb_m*periOrb_p;                     % 节点数
% ----------------------- yy_part1 (x) ------------------------
% 构造优化变量
yy_part1 = zeros(6 * n ,1) ;
for iLoop = 1 : n
    % J2000状态（未归一化）
     rvEme = BRC2Eme(aux.dep_t_jdtdb + tt_interp(iLoop)*aux.TU/86400, xx_interp(iLoop , :)', aux);
%      rvEme = crtbpRotToEme(aux.dep_t_jdtdb + tt_interp(iLoop) * aux.TU / 86400 , 3 , 10 , 'P1' , xx_interp(iLoop , :)' , aux);
    % J2000状态（归一化）
     yy_part1(6 * (iLoop - 1) + 1 : 6 * iLoop) = [rvEme(1:3) / aux.LU
       rvEme(4:6) / aux.VU];
end

% yy_part1 = xx_interp(:);
yy_part2 = diff(tt_interp);
yy_part3 = tt_interp;

aux.yFix = xx_interp(1,2);

yy = [yy_part1;
    yy_part2;
    yy_part3];

aux.node_n = n;
aux.Orbit_P = TT*periOrb_p*periOrb_m;               % 轨道总周期数


end

