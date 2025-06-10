% 分析角度对应误差
set(0,'defaultAxesFontName', 'TimesSimSun','defaultTextFontName', 'TimesSimSun');
set(0,'defaultAxesFontSize',15,'defaultTextFontSize',15)
set(0,'defaultLineLineWidth',1.5)
addpath('../subFunc(CRTBP)')
load('../DRO_SSF.mat')
r2d = 180/pi;

ii_forma = 3;

t_sample = auxSFF.naturalSFF(ii_forma).t_sample./86400;
xx_MCEMR_A_STK = auxSFF.naturalSFF(ii_forma).xx_MCEMR_A;
xx_j2kLVLH_B_STK = auxSFF.naturalSFF(ii_forma).xx_j2kLVLH_REL;
theta_A_STK = mod(atan2d(-xx_MCEMR_A_STK(:,2),xx_MCEMR_A_STK(:,1)),360);
alpha_B_STK = mod(atan2d(-xx_j2kLVLH_B_STK(:,2),xx_j2kLVLH_B_STK(:,1)),360);

alpha_B_CRTBP = theta2alpha(theta_A_STK);
theta_A_CRTBP = alpha2theta(alpha_B_STK);

theta_error = theta_A_STK-theta_A_CRTBP;
theta_error = (theta_error>=180).*(theta_error-360) + (theta_error<=-180).*(theta_error+360)+...
    ((theta_error<180)&(theta_error>-180)).*theta_error;
alpha_error = alpha_B_STK-alpha_B_CRTBP;
alpha_error = (alpha_error>=180).*(alpha_error-360) + (alpha_error<=-180).*(alpha_error+360)+...
    ((alpha_error<180)&(alpha_error>-180)).*alpha_error;

figure(1)
% plot(t_sample,theta_error)
plot(theta_A_STK,theta_error,'.')
xlabel('\theta [deg]'); ylabel('\epsilon_{\theta} [deg]')
figure(2)
% plot(t_sample,alpha_error)
plot(theta_A_STK,alpha_error,'.')
xlabel('\theta [deg]'); ylabel('\epsilon_{\alpha} [deg]')