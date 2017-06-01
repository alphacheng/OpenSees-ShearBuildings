

% K0 = 1000;              % elastic stiffness
% as_Plus = 0.04;         % strain hardening ratio for positive loading direction
% as_Neg = 0.04;          % strain hardening ratio for negative loading direction
% My_Plus = 100;          % effective yield strength for positive loading direction
% My_Neg = -100;          % effective yield strength for negative loading direction (negative value)
% Lambda_S = 5;           % Cyclic deterioration parameter for strength deterioration
% Lambda_C = 10;          % Cyclic deterioration parameter for post-capping strength deterioration
% Lambda_A = 5;           % Cyclic deterioration parameter for acceleration reloading stiffness deterioration (is not a deterioration mode for a component with Bilinear hysteretic response).
% Lambda_K = 5;           % Cyclic deterioration parameter for unloading stiffness deterioration
% c_S = 1.0;              % rate of strength deterioration. The default value is 1.0.
% c_C = 1.0;              % rate of post-capping strength deterioration. The default value is 1.0.
% c_A = 1.0;              % rate of accelerated reloading deterioration. The default value is 1.0.
% c_K = 1.0;              % rate of unloading stiffness deterioration. The default value is 1.0.
% theta_p_Plus = 0.3;     % pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
% theta_p_Neg = 0.3;      % pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (positive value)
% theta_pc_Plus = 1.5;    % post-capping rotation for positive loading direction
% theta_pc_Neg = 1.5;     % post-capping rotation for negative loading direction (positive value)
% Res_Pos = 0.2;          % residual strength ratio for positive loading direction
% Res_Neg = 0.2;          % residual strength ratio for negative loading direction (positive value)
% theta_u_Plus = 1.9;     % ultimate rotation capacity for positive loading direction
% theta_u_Neg = 1.9;      % ultimate rotation capacity for negative loading direction (positive value)
% D_Plus = 1.0;           % rate of cyclic deterioration in the positive loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
% D_Neg = 1.0;            % rate of cyclic deterioration in the negative loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
% nFactor = 0;            % elastic stiffness amplification factor, mainly for use with concentrated plastic hinge elements (optional, default = 0).

K0       = 1000; % Could be determined from drift limits
as       = 0.04; %
My       = 100;  % Could be determined from strength requirements
Lambda   = 10;   %
c        = 1.0;  %
theta_p  = 0.3;  %
theta_pc = 1.5;  %
Res      = 0.2;  %
theta_u  = 1.9;  % Related to definition of collapse?
D        = 1.0;  % Doesn't need to be determined if hysteretic behavior is assumed symmetric
nFactor  = 0;    % Leave as zero? Seems to be for explicit modeling of plastic hinges

bilinMaterial = bilinearMaterialDefinition(    K0,as,My,Lambda,c,theta_p,theta_pc,Res,theta_u,D,nFactor);
peakMaterial  = peakOrientedMaterialDefinition(K0,as,My,Lambda,c,theta_p,theta_pc,Res,theta_u,D);

peakPoints  = [0 1 -1 2 -2 3 -3 4 -4 5 -5]*0.2;
matID       = 1;
rateType    = 'StrainRate';
rateValue   = peakPoints(2)/10;

bilinObj = UniaxialMaterialAnalysis(bilinMaterial);
peakObj  = UniaxialMaterialAnalysis(peakMaterial);

bilinResults = bilinObj.runAnalysis(peakPoints,rateType,rateValue);
peakResults  =  peakObj.runAnalysis(peakPoints,rateType,rateValue);

% results_pos_env = anaobj.runAnalysis([0 max(peakPoints)],rateType,rateValue);
% results_neg_env = anaobj.runAnalysis([0 min(peakPoints)],rateType,rateValue);
% results         = anaobj.runAnalysis(         peakPoints,rateType,rateValue);


%% Make Figure
figure
hold on
% plot(results_pos_env.disp,results_pos_env.force,'k--')
% plot(results_neg_env.disp,results_neg_env.force,'k--')
plot(bilinResults.disp,bilinResults.force,'r')
plot(peakResults.disp,peakResults.force,'b')

grid on
grid minor
