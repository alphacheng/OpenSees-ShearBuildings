
function materialDefinition = bilinearMaterialDefinition(matTag, S)
% BILINEARMATERIALDEFINITION Generate an OpenSees bilinear material definition
%
%   materialDefinition = BILINEARMATERIALDEFINITION(matTag,S) uses the properties
%       in the struct S to generate a modified Ibarra-Krawinkler Bilinear
%       uniaxial material with tag matTag.
%
%   The struct S must have the following components:
%
%   S.K0
%   S.as
%   S.V_y
%   S.Lambda_S
%   S.Lambda_K
%   S.c_S
%   S.c_K
%   S.defl_p
%   S.defl_pc
%   S.Res
%   S.defl_u
%   S.D
%   S.nFactor
%

K0 = S.K0;                  % elastic stiffness
as_Plus = S.as;             % strain hardening ratio for positive loading direction
as_Neg = S.as;              % strain hardening ratio for negative loading direction
My_Plus = S.V_y;            % effective yield strength for positive loading direction
My_Neg = -S.V_y;            % effective yield strength for negative loading direction (negative value)
Lambda_S = S.Lambda_S;      % Cyclic deterioration parameter for strength deterioration
Lambda_C = S.Lambda_S;      % Cyclic deterioration parameter for post-capping strength deterioration
Lambda_A = S.Lambda_K;      % Cyclic deterioration parameter for acceleration reloading stiffness deterioration (is not a deterioration mode for a component with Bilinear hysteretic response).
Lambda_K = S.Lambda_K;      % Cyclic deterioration parameter for unloading stiffness deterioration
c_S = S.c_S;                % rate of strength deterioration. The default value is 1.0.
c_C = S.c_S;                % rate of post-capping strength deterioration. The default value is 1.0.
c_A = S.c_K;                % rate of accelerated reloading deterioration. The default value is 1.0.
c_K = S.c_K;                % rate of unloading stiffness deterioration. The default value is 1.0.
theta_p_Plus = S.defl_p;    % pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
theta_p_Neg = S.defl_p;     % pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (positive value)
theta_pc_Plus = S.defl_pc;  % post-capping rotation for positive loading direction
theta_pc_Neg = S.defl_pc;   % post-capping rotation for negative loading direction (positive value)
Res_Pos = S.Res;            % residual strength ratio for positive loading direction
Res_Neg = S.Res;            % residual strength ratio for negative loading direction (positive value)
theta_u_Plus = S.defl_u;    % ultimate rotation capacity for positive loading direction
theta_u_Neg = S.defl_u;     % ultimate rotation capacity for negative loading direction (positive value)
D_Plus = S.D;               % rate of cyclic deterioration in the positive loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
D_Neg = S.D;                % rate of cyclic deterioration in the negative loading direction (this parameter is used to create asymmetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
nFactor = S.nFactor;        % elastic stiffness amplification factor, mainly for use with concentrated plastic hinge elements (optional, default = 0).

materialDefinition = sprintf('uniaxialMaterial Bilin %i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',...
    matTag,K0,as_Plus,as_Neg,My_Plus,My_Neg,Lambda_S,Lambda_C,Lambda_A,Lambda_K,...
    c_S,c_C,c_A,c_K,theta_p_Plus,theta_p_Neg,theta_pc_Plus,theta_pc_Neg,...
    Res_Pos,Res_Neg,theta_u_Plus,theta_u_Neg,D_Plus,D_Neg,nFactor);

end
