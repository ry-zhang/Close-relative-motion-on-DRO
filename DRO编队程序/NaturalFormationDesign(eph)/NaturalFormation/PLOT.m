
load("r_MCRLVLH_rel_ref.mat");
figure(1)
plot3(r_MCRLVLH_rel_ref(1,1),r_MCRLVLH_rel_ref(1,2),r_MCRLVLH_rel_ref(1,3),'g^');hold on;
plot3(r_MCRLVLH_rel_ref(end,1),r_MCRLVLH_rel_ref(end,2),r_MCRLVLH_rel_ref(end,3),'r*');hold on;
plot3(r_MCRLVLH_rel_ref(:,1),r_MCRLVLH_rel_ref(:,2),r_MCRLVLH_rel_ref(:,3),'b');hold on;