function aux = SaveOrbit2Aux(aux,satA,satB,t0,tf,tstep)
% 将A星与B星的轨道数据保存至结构体

% A星
xx_j2k_A = PosVel_j2k(satA,t0,tf,tstep);
xx_MCEMR_A = PosVel_MCEMR(satA,t0,tf,tstep);
aux.x0_MCEMR_A = xx_MCEMR_A(1,2:end);
aux.x0_j2k_A = xx_j2k_A(1,2:end);
aux.t_sample = xx_MCEMR_A(:,1);
aux.xx_MCEMR_A = xx_MCEMR_A(:,2:end);
aux.xx_j2k_A = xx_j2k_A(:,2:end);
% B星
xx_j2k_B = PosVel_j2k(satB,t0,tf,tstep);
xx_MCEMR_B = PosVel_MCEMR(satB,t0,tf,tstep);
aux.x0_j2k_B = xx_j2k_B(1,2:end);
aux.x0_MCEMR_B = xx_MCEMR_B(1,2:end);
aux.xx_j2k_B = xx_j2k_B(:,2:end);
aux.xx_MCEMR_B = xx_MCEMR_B(:,2:end);
% 相对运动(B星-A星)
xx_j2kLVLH_REL = PosVel_j2kLVLH(satB,satA,t0,tf,tstep);
xx_MCEMRLVLH_REL = PosVel_MCEMRLVLH(satB,satA,t0,tf,tstep);
aux.x0_j2kLVLH_REL = xx_j2kLVLH_REL(1,:);
aux.x0_MCEMRLVLH_REL = xx_MCEMRLVLH_REL(1,:);
aux.xx_j2kLVLH_REL = xx_j2kLVLH_REL;
aux.xx_MCEMRLVLH_REL = xx_MCEMRLVLH_REL;