function psi21 = zsingularTM(dt,t0,para,coe)
% 状态转移矩阵求z方向转移奇点

tf = t0+dt;

theta1_0 = 0; theta2_0 = 0;
k0 = 0; k1 = 0; k2 = 1; 
rel_motion_temp_t01 = generalSol_relMotion(t0,k0,k1,k2,theta1_0,theta2_0,para,coe);
rel_motion_temp_t02 = generalSol_relMotion(t0,k0,k1,k2,theta1_0,theta2_0+pi/2,para,coe);
Phi_t0 = [rel_motion_temp_t01([3,6]),rel_motion_temp_t02([3,6])];

rel_motion_temp_tf1 = generalSol_relMotion(tf,k0,k1,k2,theta1_0,theta2_0,para,coe);
rel_motion_temp_tf2 = generalSol_relMotion(tf,k0,k1,k2,theta1_0,theta2_0+pi/2,para,coe);
Phi_tf = [rel_motion_temp_tf1([3,6]),rel_motion_temp_tf2([3,6])];
Psi = Phi_tf*Phi_t0^(-1);

psi21 = Psi(1,2);


