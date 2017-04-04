

K0 = 1000;              % elastic stiffness
as_Plus = 0.04;         % strain hardening ratio for positive loading direction
as_Neg = 0.04;          % strain hardening ratio for negative loading direction
My_Plus = 100;          % effective yield strength for positive loading direction
My_Neg = -100;          % effective yield strength for negative loading direction (negative value)
Lambda_S = 5;           % Cyclic deterioration parameter for strength deterioration
Lambda_C = 5;          % Cyclic deterioration parameter for post-capping strength deterioration
Lambda_A = 5;           % Cyclic deterioration parameter for acceleration reloading stiffness deterioration (is not a deterioration mode for a component with Bilinear hysteretic response).
Lambda_K = 5;           % Cyclic deterioration parameter for unloading stiffness deterioration
c_S = 1.0;              % rate of strength deterioration. The default value is 1.0.
c_C = 1.0;              % rate of post-capping strength deterioration. The default value is 1.0.
c_A = 1.0;              % rate of accelerated reloading deterioration. The default value is 1.0.
c_K = 1.0;              % rate of unloading stiffness deterioration. The default value is 1.0.
theta_p_Plus = 0.3;     % pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
theta_p_Neg = 0.3;      % pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (positive value)
theta_pc_Plus = 1.5;    % post-capping rotation for positive loading direction
theta_pc_Neg = 1.5;     % post-capping rotation for negative loading direction (positive value)
Res_Pos = 0.2;          % residual strength ratio for positive loading direction
Res_Neg = 0.2;          % residual strength ratio for negative loading direction (positive value)
theta_u_Plus = 1.9;     % ultimate rotation capacity for positive loading direction
theta_u_Neg = 1.9;      % ultimate rotation capacity for negative loading direction (positive value)
D_Plus = 1.0;           % rate of cyclic deterioration in the positive loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
D_Neg = 1.0;            % rate of cyclic deterioration in the negative loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.



materialDefinition = sprintf('uniaxialMaterial ModIMKPeakOriented 1 %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',...
    K0,as_Plus,as_Neg,My_Plus,My_Neg,Lambda_S,Lambda_C,Lambda_A,Lambda_K,...
    c_S,c_C,c_A,c_K,theta_p_Plus,theta_p_Neg,theta_pc_Plus,theta_pc_Neg,...
    Res_Pos,Res_Neg,theta_u_Plus,theta_u_Neg,D_Plus,D_Neg);

peakPoints  = [0 1 -1 2 -2 3 -3 4 -4 5 -5]*0.2;
matID       = 1;
rateType    = 'StrainRate';
rateValue   = peakPoints(2)/10;


anaobj = UniaxialMaterialAnalysis(materialDefinition);
%anaobj.deleteFilesAfterAnalysis = false;
results_pos_env = anaobj.runAnalysis([0 max(peakPoints)],rateType,rateValue);
results_neg_env = anaobj.runAnalysis([0 min(peakPoints)],rateType,rateValue);
results         = anaobj.runAnalysis(         peakPoints,rateType,rateValue);


%% Make Figure
figure
hold on
plot(results_pos_env.disp,results_pos_env.force,'k--')
plot(results_neg_env.disp,results_neg_env.force,'k--')
plot(results.disp,results.force,'r')

grid on
grid minor
