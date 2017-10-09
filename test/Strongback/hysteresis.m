
clear all; close all; home; %#ok<CLALL>

% --- Springs ---
spring.K0 = 1807.3;
spring.as       =  0.03;  % strain hardening ratio
spring.V_y = 343.38;
spring.Lambda_S =  200.00;  % Cyclic deterioration parameter - strength
spring.Lambda_K =  0.00;  % Cyclic deterioration parameter - stiffness
spring.c_S      =  1.00;  % rate of deterioration - strength
spring.c_K      =  1.00;  % rate of deterioration - stiffness
spring.defl_p   =  9e99;
spring.defl_pc  =  9e99;
spring.defl_u   =  9e99;
spring.Res      =  0.30;  % residual strength ratio (relative to yield)
spring.D        =  1.00;  % rate of cyclic deterioration
spring.nFactor  =  0.00;  % elastic stiffness amplification factor

spring.definition = bilinearMaterialDefinition(1, spring);



%================================== Run stuff =================================%
anaobj = UniaxialMaterialAnalysis(spring.definition);

peakPoints  = [0 11 -11 12 -12 13 -13 14 -14 15 -15]/5;
rateType    = 'StrainRate';
rateValue   = peakPoints(2)/10;

hysteretic_pos_env = anaobj.runAnalysis([0 1.2*max(peakPoints)],rateType,rateValue);
hysteretic_neg_env = anaobj.runAnalysis([0 1.2*min(peakPoints)],rateType,rateValue);
hysteretic         = anaobj.runAnalysis(           peakPoints,  rateType,rateValue);

figure
hold on
plot(hysteretic_pos_env.disp, hysteretic_pos_env.force, 'k--')
plot(hysteretic_neg_env.disp, hysteretic_neg_env.force, 'k--')
plot(hysteretic.disp, hysteretic.force, 'r')
grid on
title('Hysteretic behavior of spring')
xlabel('Deflection (ft)')
ylabel('Force (kip)')
