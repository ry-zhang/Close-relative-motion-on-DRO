%% natural and forced 2-sat formations in DRO based on linear relative motion of CR3BP

%% natural formation
clear
load('NaturalFormationData.mat')
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);

% 标称轨道及相对运动 的积分
dt_natural = 2*para.T0; % 积分时间
sol_natural = ode113(@(t,x)eom_rel3b(t,x,para.mu),[0 dt_natural], [x0_DRO, x0_REL], opts_ode);
abs_motion_natural = sol_natural.y(1:6,:);
Natural_Traj = sol_natural.y(7:12,:);
Natural_Traj_km = Natural_Traj*para.r_norma;

%% plot (natural formation)
figure(314)
plot3(Natural_Traj_km(1,1),Natural_Traj_km(2,1),Natural_Traj_km(3,1),'g^'); hold on
plot3(Natural_Traj_km(1,:),Natural_Traj_km(2,:),Natural_Traj_km(3,:),'Color',[237, 177, 32]/255,'LineWidth',1.5);
plot3(Natural_Traj_km(1,end),Natural_Traj_km(2,end),Natural_Traj_km(2,end),'rv'); 
plot3(0,0,0,'ks'); 
text(2,0.8,0,'\leftarrow Leader','FontSize',13)
hold off
xlabel('x_L [km]'); ylabel('y_L [km]'); zlabel('z_L [km]');
legend('Initial Position','Final Position','Natural Trajectory')
axis equal; title('Natural Formation'); set(gca,'FontSize',13)
grid on; grid minor
box on
view(0,90)
xlim([-20,20])
ylim([-60,10])

%% forced formation
clear
load('ForcedFormationData')
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-20);

Forced_Traj_all = [];
num_dv = size(dv_all,1);

for ii_index = 1:num_dv
%     disp(ii_index)
    x0_DRO_loop = x0_DRO_all(ii_index,:);
    x0_REL_loop = x0_REL_all(ii_index,:);
    [~,y_temp] = ode113(@(t,x)eom_rel3b(t,x,para.mu),[0 dt], [x0_DRO_loop, x0_REL_loop], opts_ode);
    Forced_Traj_all = [Forced_Traj_all;y_temp];
end

r0_REL_all_km = x0_REL_all(:,1:3)*para.r_norma;
r0_DRO_all_km = x0_DRO_all(:,1:3)*para.r_norma;
Forced_Traj_all_km = Forced_Traj_all*para.r_norma;

%% plot (forced formation)
fig = figure(315);
hold off
plot3(r0_REL_all_km(1,1),r0_REL_all_km(1,2),r0_REL_all_km(1,3),'g^');  hold on
plot3(r0_REL_all_km(2:num_dv,1),r0_REL_all_km(2:num_dv,2),r0_REL_all_km(2:num_dv,3),'bo');
plot3(r0_REL_all_km(end,1),r0_REL_all_km(end,2),r0_REL_all_km(end,3),'rv');
plot3(Forced_Traj_all_km(:,7),Forced_Traj_all_km(:,8),Forced_Traj_all_km(:,9),'Color',[0, 114, 189]/255,'LineWidth',1.5); 
plot3(0,0,0,'ks'); 
text(0.05,0.03,0,'\leftarrow Leader','FontSize',13)
hold off
xlabel('x_L [km]'); ylabel('y_L [km]'); zlabel('z_L [km]');
legend('Initial Position','Patch Position','Final Position','Forced Trajectory')
axis equal; title('Forced Formation'); set(gca,'FontSize',13)
grid on; grid minor
box on
view(0,90)
xlim([-1.5,1.5])
ylim([-1.5,1.5])