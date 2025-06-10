function [yy,aux] = getIC_periOrb(aux)
%
% ������ڹ����ֵ
%
% ���ߣ�wm
% 2021��11��3��
%%%%%%%%%%%%%%%%%%%%%%%%
% ������
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

% ���������ڷ�Χ
Pmin = min(X0Mtx(: , end - 1));
Pmax = max(X0Mtx(: , end - 1));
% ������������Χ
Cmin = min(X0Mtx(: , end));
Cmax = max(X0Mtx(: , end));
fprintf('��������ͣ�%s \n' , aux.periOrb_Modi)
fprintf('��������ڣ�[%0.4f , %0.4f] \n' , Pmin , Pmax )
fprintf('�����������[%0.4f , %0.4f] \n' , Cmin , Cmax)
fprintf('--------------------------------------------------- \n');
% �������ڱ���
periOrb_m = aux.periOrb_m;

% ���ڹ��һȦ�ڵ���
periOrb_n = aux.periOrb_n;
% �ж��Ƿ�����
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
        fprintf('��ǰ���ڶ�Ӧ���������! \n');
        fprintf('�������������ڷ�Χ��[%0.4f , %0.4f],' , Pmin, Pmax )
        prompt = ('���������������ڣ� \n');
        num = input(prompt);
        if num<=Pmax && num>=Pmin
            TT = num;
        else
            fprintf('������󣬳�����ֹ! \n');
            return;
        end
    end
    % �Բ�ֵ���д��У��
    fprintf(' ------- �����ֵУ�� --------- \n')
    x0 =  interp1(Orbit_data(:,7) , Orbit_data(:,1:6) , TT , 'spline');  %    �����ʼλ��
    x0 = Refine_x0(x0, TT, aux);
    options = odeset('RelTol',1e-13,'AbsTol',1e-16);
    [tt, xx]  = ode113(@crtbpEqm3D, [0, TT*periOrb_m], x0 , options, aux);
    figure (100)
    plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'b');
    % ����������ֵ
    tt_interp = linspace(0 , TT*periOrb_m , periOrb_n*periOrb_m*periOrb_p)';
    xx_interp = interp1(tt , xx , tt_interp' , 'spline') ;
%     fprintf('��ʼ̫����ǣ�%0.0f�� \n' , rad2deg(aux.Sun_angle))
    fprintf('���ת����%0.0f����ʼ���� %0.4f\n' , periOrb_m*periOrb_p,TT*periOrb_m)
else
   
    
    if isempty(aux.periOrb_T)
        % �ж� p,q�Ƿ�Լ
        [~, gcd] = is_coprime(aux.resoOrb_p  , aux.resoOrb_q);
        
        % ����ȣ�p : q
        periOrb_p = aux.periOrb_p/gcd ;
        periOrb_q  = aux.periOrb_q/gcd ;
        TT = 2*pi*periOrb_q/periOrb_p;                             % ��Ȧ�������
        fprintf('%0.0f ��%0.0f ��Ȧ������ڣ�%0.4f \n' , periOrb_p , periOrb_q,TT)
%         TT = 1/abs(aux.oms)*2*pi*periOrb_q/periOrb_p;
%         % �жϸù���ȹ���Ƿ����
%         fprintf('����̫�����ڹ����%0.0f ��%0.0f ��Ȧ������ڣ�%0.4f \n' , periOrb_p , periOrb_q,TT)
    else
        periOrb_p = 1;
        %         periOrb_q = 1;
        TT = str2double(aux.periOrb_T);
        fprintf('%0.0f Ȧ������ڣ�%0.4f \n' , periOrb_p ,TT)
        %      TT = 1/abs(aux.oms)*TT;
        %     fprintf('����������̫�����ڹ����%0.0f ��%0.0f ��Ȧ�����������Ϊ��%0.4f \n' , periOrb_p , periOrb_q,TT)
        
    end
    
    
    % ��ʼʱ��̫�����
    
    
    % ƥ�䵥Ȧ����
    
    % ע�⣬������˷���չ���
    
    
    if TT>max(X0Mtx(:,7)) || TT<min(X0Mtx(:,7))
        periOrb_p = 1;                                              % ���漰�����
%         periOrb_q = 1;
        fprintf('--------------------------------------------------- \n');
        fprintf('��ǰ���ڶ�Ӧ���������! \n');
        fprintf('�������������ڷ�Χ��[%0.4f , %0.4f],' , Pmin, Pmax )
        prompt = ('���������������ڣ� \n');
        num = input(prompt);
        if num<=Pmax && num>=Pmin
            TT = num;
        else
            fprintf('������󣬳�����ֹ! \n');
            return;
        end
    end
%     Sun_angle0 = aux.Sun_angle;
%     fprintf('��ʼ̫����ǣ�%0.0f�� \n' , rad2deg(Sun_angle0))
    fprintf('���ת����%0.0f����ʼ���� %0.4f\n' , periOrb_m*periOrb_p,TT*periOrb_p*periOrb_m)
    Orbit_data = X0Mtx;
    % % �жϴ������������
    % if strcmp(aux.orb_type , 'DRO')
    %     Orbit_data = aux.DRO;
    % end
    
    x0 =  interp1(Orbit_data(:,7) , Orbit_data(:,1:6) , TT , 'spline');  %    �����ʼλ��
  
    % �Բ�ֵ���д��У��
    fprintf(' ------- �����ֵУ�� --------- \n')
    x0 = Refine_x0(x0, TT, aux);
    Tspan = [0, TT*periOrb_p*periOrb_m];
    options = odeset('RelTol',1e-13,'AbsTol',1e-16);
    [tt, xx]  = ode113(@crtbpEqm3D, Tspan, x0 , options, aux);
    
    % ����������ֵ
    tt_interp = linspace(0 , TT*periOrb_p*periOrb_m , periOrb_n*periOrb_m*periOrb_p)';
    xx_interp = interp1(tt , xx , tt_interp' , 'spline') ;
    
    
    
    % �����Ż�����
end
% �����ԡ�
figure(1); hold on; grid on;
crtbpMarkEM;
plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'b');
plot3(xx_interp(: , 1) , xx_interp(: , 2) , xx_interp(: , 3) , 'bo');
    
% xx_interp = xx_interp';
n = periOrb_n*periOrb_m*periOrb_p;                     % �ڵ���
% ----------------------- yy_part1 (x) ------------------------
% �����Ż�����
yy_part1 = zeros(6 * n ,1) ;
for iLoop = 1 : n
    % J2000״̬��δ��һ����
     rvEme = BRC2Eme(aux.dep_t_jdtdb + tt_interp(iLoop)*aux.TU/86400, xx_interp(iLoop , :)', aux);
%      rvEme = crtbpRotToEme(aux.dep_t_jdtdb + tt_interp(iLoop) * aux.TU / 86400 , 3 , 10 , 'P1' , xx_interp(iLoop , :)' , aux);
    % J2000״̬����һ����
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
aux.Orbit_P = TT*periOrb_p*periOrb_m;               % �����������


end

