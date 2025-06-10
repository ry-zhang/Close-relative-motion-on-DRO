%% 计算DRO非线性相对运动模型的非对称性
%2023-12-11
% by ZhangRuYue
% email:zhangruyue22@csu.ac.cn
close all
clear
addpath('../subF_eom(CR3BP)')

%% 常数与变量
mu_E = 398600.44; % km^3*s^-2
mu_M = 4904.8695; % km^3*s^-2  //4902.56146783419
% con.mu = 0.01211; % 
con.mu = 0.01215; % 20200531
con.r_norma = 3.84399*10^5; % km
%T_M = 27.321661; % day
con.T_norma =2*pi*sqrt(con.r_norma^3/(mu_E+mu_M)); % s
con.T_norma_day = con.T_norma/3600/24;
con.v_norma = con.r_norma/con.T_norma; % km/s
%con.v_norma = con.r_norma/(T_M*24*2600); % km/s
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);
T_DRO=1/2*con.T_norma_day;
isplot = 1;

%% 积分与计算
% 绝对运动初值
% 初值由共同质心旋转坐标系S下的动力学（陈冠华程序）导出：原点位于共同质心，x轴从地球指向月球，y轴指向月球绕地球旋转方向
% 绝对动力学位于月球旋转坐标系M：原点位于月球，x轴从地球指向月球，y轴指向月球绕地球旋转的方向
% 惯性坐标系：原点位于共同质心，x、y、z轴与初始时刻的M坐标系相同
% 从S至M的位置转换：rM = rS + [1-mu;0]; 速度转换： vM = vS;
% load('InitialState_DRO_2to1')
% load('DRO_all.mat')
% load('DRO_all_SEsyPeriod')
load('DRO_all_MoonPeriod')
%  A=[4];
A=[ 1, 4, 6,7 ,8]
for jj_index =2:2
 numDRO = A(1,jj_index);%4 6 7 8 9 /2:1 3:1 4:1 5:1 6:1
% numDRO = jj_index;%4 6 7 8 9 /2 3 4 5 6
x0_DRO_S_2d = state_ini_all(numDRO,:);
x0_DRO_M_3d = [(x0_DRO_S_2d(1:2)-[1-con.mu,0]),0,x0_DRO_S_2d(3:4),0];
%x0_DRO_S_3d = [x0_DRO_S_2d(1:2),0,x0_DRO_S_2d(3:4),0];
para.T0 = J_period_all(numDRO,2)*2*pi; % 轨道周期
para.T01= J_period_all(numDRO,2)*27.284429211881122*24*60*60;
%% %% 线性化相对运动的周期状态转移矩阵（φ(0)^-1 φ(T)）
sol_temp = ode113(@(t,x)eomM_rel3b(t,x,con.mu),[0 1*para.T0], [x0_DRO_M_3d, zeros(1,6), reshape(eye(6),1,36)], opts_ode);
M_REL_lin = reshape(sol_temp.y(13:end,end),6,6);
[Meigve, Meigva] = eig(M_REL_lin);
% remove small numbers in Meigve and obtain the accurate eigen vector
Meigve_imga = imag(Meigve); Meigve_imga(abs(Meigve_imga)<1e-5) = 0;
Meigve_real = real(Meigve); Meigve_real(abs(Meigve_real)<1e-5) = 0;
Meigve = Meigve_real+Meigve_imga*1i;

% 判断单位特征值在第1-2行，还是第3-4行
[a,sort_index] = sort(real(diag(Meigva(1:4,1:4))));
index_temp1 = 3*(sort_index(1)>2) + 1*(sort_index(1)<=2);
index_temp2 = 1*(sort_index(1)>2) + 3*(sort_index(1)<=2);
H = blkdiag([real(Meigva(index_temp1,index_temp1)),-imag(Meigva(index_temp1,index_temp1)); ...
    imag(Meigva(index_temp1,index_temp1)),real(Meigva(index_temp1,index_temp1))],...
    [1,1; 0,1],...
    [real(Meigva(5,5)),-imag(Meigva(5,5)); imag(Meigva(5,5)),real(Meigva(5,5))]);
p2_real = real(Meigve(:,index_temp1)); p2_imag = imag(Meigve(:,index_temp1+1));
p3 = real(Meigve(:,index_temp2)); p3 = sign(p3(2)).*p3;
p4 = pinv(M_REL_lin-eye(6))*p3; % 由于M_REL_lin-eye(6)奇异，采用广义逆解
try
    p4 = p4-p4(2)/p3(2)*p3;
catch
end
p3 = p3/norm(p3); 
if norm(p4)>10 % 精度太高有时候算出来的p4超级大
    p4 = pinv(M_REL_lin-eye(6),1e-10)*p3;
end
p4 = p4/norm(p3);

p6_real = real(Meigve(:,6)); p6_imag = imag(Meigve(:,6));

S = [p2_real,p2_imag,p3,p4,p6_real,p6_imag];
flag_Jordan = norm(S*H*S^-1-M_REL_lin);

J = 1/para.T0*logm(H);
J = J.*(abs(J)>1e-8);
omega1 = J(2,1);
omega2 = J(3,4); % 正好是1/T0
omega3 = J(6,5);

expJt = @(t)blkdiag([cos(omega1*t),-sin(omega1*t);sin(omega1*t),cos(omega1*t)],...
    [1,omega2*t;0,1],...
    [cos(omega3*t),-sin(omega3*t);sin(omega3*t),cos(omega3*t)]);
flag_Floquet = norm(M_REL_lin * S * expJt(-para.T0) - S);

% resort data
Meigva_diag = diag(Meigva)';
Meigva_diag = Meigva_diag([index_temp1,index_temp1+1,index_temp2,index_temp2+1,5,6]);
Meigve = Meigve(:,[index_temp1,index_temp1+1,index_temp2,index_temp2+1,5,6]);

% alpha其实就是算拟周期的时候的旋转量，计算的时候取正虚部的特征值计算
alpha1 = [atan2(imag(Meigva_diag(2)),real(Meigva_diag(2))),...
        atan2(imag(Meigva_diag(6)),real(Meigva_diag(6)))];
T12_propotion = 2*pi./alpha1;
% save data

Sol_linear.vec1 = p2_real;
Sol_linear.vec2 = p2_imag;
Sol_linear.vec3 = p3;
Sol_linear.vec4 = p4;
Sol_linear.vec5 = p6_real;
Sol_linear.vec6 = p6_imag;

% 这里采用第二个和第六个特征值，因为其虚部是正的，如此可以确保alpha大于0
EigenVector.p1_real = p2_real;
EigenVector.p1_imag = p2_imag;
EigenVector.p3 = p3;
EigenVector.p4 = p4;
EigenVector.p5_real = p6_real;
EigenVector.p5_imag = p6_imag;
save ../subF_eom(CR3BP)/FloquetEig12 Meigve Meigva_diag EigenVector Sol_linear state_ini_all J_period_all x0_DRO_M_3d para con


%% 非对称度图
% if isplot == 1
    dt =1*para.T0;
    length_t = 2000;
%     t_sample = linspace(0,dt,length_t);
    t_sample = linspace(0,dt,length_t);
%     t_sample=1:60:dt;
    opts = odeset('RelTol',1e-13,'AbsTol',1e-20); 
%      opts = odeset('RelTol',1e-1,'AbsTol',1e-1); 
        eval(['x0_REL = Sol_linear.vec',num2str(3),'*1e-5;'])
        m(jj_index)=50/x0_REL(2)/con.r_norma;%前面的数字代表真实编队尺度km
        x0_REL= m(jj_index)* x0_REL';
        x0_DRO_M_3d1=[x0_DRO_M_3d(1:3)*con.r_norma,x0_DRO_M_3d(4:6)*con.v_norma];
        x0_REL1=[x0_REL(1:3)*con.r_norma,x0_REL(4:6)*con.v_norma];
%         mu_con=4.0351e+05*con.mu;

        % 标称轨道及相对运动 的积分
        sol = ode113(@(t,x)eom_rel3b_alpha(t,x,con.mu),[0 dt], [x0_DRO_M_3d,x0_REL], opts);
% %      sol = ode113(@(t,x)eom_rel3b_alphacon(t,x,con.mu),[0 dt], [x0_DRO_M_3d1,x0_REL1], opts);
        sol_sample = deval(sol,t_sample);
        abs_motion_M= sol_sample(1:6,:);
        abs_motion_M_con =[sol_sample(1:3,:)*con.r_norma;sol_sample(4:6,:)*con.v_norma];
        rel_motion_L = [sol_sample(7:9,:);sol_sample(10:12,:)];
        rel_motion_L_con= [sol_sample(7:9,:)*con.r_norma;sol_sample(10:12,:)*con.v_norma];%L
%         abs_pos=sqrt(abs_motion_M(1,:).^2+abs_motion_M(2,:).^2+abs_motion_M(3,:).^2);
%         abs_vel=sqrt(abs_motion_M(4,:).^2+abs_motion_M(5,:).^2+abs_motion_M(6,:).^2);
%         rel_pos=sqrt(rel_motion_L(1,:).^2+rel_motion_L(2,:).^2+rel_motion_L(3,:).^2);
%         rel_vel=sqrt(rel_motion_L(4,:).^2+rel_motion_L(5,:).^2+rel_motion_L(6,:).^2);
        abs_pos=sqrt(abs_motion_M_con(1,:).^2+abs_motion_M_con(2,:).^2+abs_motion_M_con(3,:).^2);
        abs_vel=sqrt(abs_motion_M_con(4,:).^2+abs_motion_M_con(5,:).^2+abs_motion_M_con(6,:).^2);
        rel_pos=sqrt(rel_motion_L_con(1,:).^2+rel_motion_L_con(2,:).^2+rel_motion_L_con(3,:).^2);
        rel_vel=sqrt(rel_motion_L_con(4,:).^2+rel_motion_L_con(5,:).^2+rel_motion_L_con(6,:).^2);

        rel_motion_L1=rel_motion_L_con;
%         rel_motion_L1=rel_motion_L;
%         X_Xr(1:6,:)=abs_motion_M;
%         X_Xr(7:12,:)=rel_motion_L;
        X_Xr(1:6,:)=abs_motion_M_con;
        X_Xr(7:12,:)=rel_motion_L_con;
         
        % 将绝对运动轨道转移至惯性系中 m
        
         for ii = 1:2000
               Al(1:4,ii)=alphacon(X_Xr(1:3,ii),X_Xr(4:6,ii),X_Xr(7:9,ii),X_Xr(10:12,ii));
%             Al(jj_index,ii)=alphacon(X_Xr(1:3,ii),X_Xr(4:6,ii),X_Xr(7:9,ii),X_Xr(10:12,ii));
         end
 %-------------------修改-------------
        rel_motion_L1=T_TCO2TCR_CR3BP(rel_motion_L_con',abs_motion_M_con','LVLH',con.mu)';%M
%          注销输出L系，不注销M系
 %-------------------修改-------------

%          rel_motion_L1= T_TCO2TCR_CR3BP(rel_motion_L',abs_motion_M','LVLH',con.mu)';%M

       
%         idx=Al(1,:)<0.4;
%         s=Al(1,:);
       
%          rel_motion_M_linear = T_TCR2TCO_CR3BP(rel_motion_M_linear',abs_motion_M','VVLH',con.mu)';

%         rel_motion_L_linear =T_TCO2TCR_CR3BP(rel_motion_L_linear',abs_motion_M','LVLH',con.mu)';
%         rel_motion_VVLH_linear = T_TCR2TCO_CR3BP(rel_motion_M_linear',abs_motion_M','VVLH',con.mu)';

       figure(1)
            
            plot(rel_motion_L1(1,1),rel_motion_L1(2,1),'g^'); hold on
            plot(rel_motion_L1(1,end),rel_motion_L1(2,end),'r*');  hold on

%              plot(rel_motion_L_Nolinear(1,:),rel_motion_L_Nolinear(2,:),'Color',[0, 114, 189]/255,'LineWidth',1.5); hold on
            
%           Al=alpha(abs_motion_L1_M',rel_motion_L1_linear',con.mu);

% %         else
%             plot(t_sample/para.T0,rel_motion_L_Nolinear(3,:),'LineWidth',1.5)
%             xlim([min(t_sample/para.T0), max(t_sample/para.T0)]) 


%         x=rel_motion_L1(1,2:1999);
%         y=rel_motion_L1(2,2:1999);
%         z=Al(jj_index,2:1999);
        x=rel_motion_L1(1,:);
        y=rel_motion_L1(2,:);
        z=Al(1,:)*100;
        %z=s(idx);
        scatter(x,y,8,z,'filled'); hold on
       
end
        legend('初始时刻','终端时刻')
        xlabel('\itx \rm[km]'); ylabel('\ity \rm[km]')
       
        axis equal;
        title('月心旋转系下 50km编队尺度 相对运动非对称度'); 
        grid on; grid minor
       
        h=colorbar('Direction','reverse');
        set(get(h,'Title'),'string','%');
   
  xlist=linspace(1,14,14);
  figure(2)
  subplot(3,2,1)
  plot(t_sample/para.T0,rel_pos); hold on
  title('星间相对距离');    
  xlabel('\it相位 \rm[T]'); ylabel('\itr \rm[km]')
  subplot(3,2,2)
  plot(t_sample/para.T0,rel_vel*1000); hold on
  title('星间相对速度');  
  xlabel('\it相位 \rm[T]'); ylabel('\itv \rm[m/s]')
  subplot(3,2,3)
  plot(t_sample/para.T0,abs_pos/10000); hold on
  title('abspos');    
  xlabel('\it相位 \rm[T]'); ylabel('\itr \rm[万km]')%距离月心
  subplot(3,2,4)
  plot(t_sample/para.T0,abs_vel*1000); hold on
  title('absvel');  
%   xlabel('\it相位 \rm[T]'); ylabel('\itv \rm[m/s]')
%   subplot(2,2,3)
%   plot(t_sample/para.T0,Al(1,:)*100) ;hold on
%   title('非对称度');
%   xlabel('\it相位 \rm[T]'); ylabel('\it \rm[%]')
%   subplot(3,2,3)
%   plot(t_sample/para.T0,Al(1,:)) ;hold on
%   xlabel('\it相位 \rm[T]'); ylabel('\ita \rm[m/s2]')
%   title('非对称度');
%   subplot(3,2,4)
%   plot(t_sample/para.T0,Al(2,:)*1000) ;hold on
%   xlabel('\it相位 \rm[T]'); ylabel('\ita \rm[m/s2]')
%   title('非对称相对加速度');
  subplot(3,2,5)
  plot(t_sample/para.T0*T_DRO,Al(2,:)) ;hold on
  xlabel('\it时间 \rm[T]'); ylabel('\ita \rm[m/s2]')
%  xticks(xlist)
%   xticklabels(string(xlist))
  title('对称相对加速度');
  subplot(3,2,6)
  plot(t_sample/para.T0*T_DRO,Al(1,:)*100) ;hold on
  xlabel('\it时间 \rm[T]'); ylabel('\ita \rm[%]')
%   xticks(xlist)
%   xticklabels(string(xlist))
  title('非对称度');
%   subplot(4,2,8)
%   plot(t_sample/para.T0,Al(4,:)*1000) ;hold on
%   xlabel('\it相位 \rm[T]'); ylabel('\ita \rm[m/s2]')
%   title('|a21|+|a3|');
%      
    
 
%   hold off



