
function materialDefinition = bilinearMaterialDefinition(matTag, K0, as, My, Lambda, c, theta_p, theta_pc, Res, theta_u, D, nFactor)

%K0;                        % elastic stiffness
as_Plus = as;               % strain hardening ratio for positive loading direction
as_Neg = as;                % strain hardening ratio for negative loading direction
My_Plus = My;               % effective yield strength for positive loading direction
My_Neg = -My;               % effective yield strength for negative loading direction (negative value)
Lambda_S = Lambda;          % Cyclic deterioration parameter for strength deterioration
Lambda_C = Lambda;          % Cyclic deterioration parameter for post-capping strength deterioration
Lambda_A = Lambda;          % Cyclic deterioration parameter for acceleration reloading stiffness deterioration (is not a deterioration mode for a component with Bilinear hysteretic response).
Lambda_K = Lambda;          % Cyclic deterioration parameter for unloading stiffness deterioration
c_S = c;                    % rate of strength deterioration. The default value is 1.0.
c_C = c;                    % rate of post-capping strength deterioration. The default value is 1.0.
c_A = c;                    % rate of accelerated reloading deterioration. The default value is 1.0.
c_K = c;                    % rate of unloading stiffness deterioration. The default value is 1.0.
theta_p_Plus = theta_p;     % pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
theta_p_Neg = theta_p;      % pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (positive value)
theta_pc_Plus = theta_pc;   % post-capping rotation for positive loading direction
theta_pc_Neg = theta_pc;    % post-capping rotation for negative loading direction (positive value)
Res_Pos = Res;              % residual strength ratio for positive loading direction
Res_Neg = Res;              % residual strength ratio for negative loading direction (positive value)
theta_u_Plus = theta_u;     % ultimate rotation capacity for positive loading direction
theta_u_Neg = theta_u;      % ultimate rotation capacity for negative loading direction (positive value)
D_Plus = D;                 % rate of cyclic deterioration in the positive loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
D_Neg = D;                  % rate of cyclic deterioration in the negative loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
%nFactor;                   % elastic stiffness amplification factor, mainly for use with concentrated plastic hinge elements (optional, default = 0).

materialDefinition = sprintf('uniaxialMaterial Bilin %i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',...
    matTag,K0,as_Plus,as_Neg,My_Plus,My_Neg,Lambda_S,Lambda_C,Lambda_A,Lambda_K,...
    c_S,c_C,c_A,c_K,theta_p_Plus,theta_p_Neg,theta_pc_Plus,theta_pc_Neg,...
    Res_Pos,Res_Neg,theta_u_Plus,theta_u_Neg,D_Plus,D_Neg,nFactor);

end
