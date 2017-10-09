classdef DeterioratingSpring
% DETERIORATINGSPRING Class for OpenSees modified Ibarra-Medina-Krawinkler materials

properties
    matTag          % Material tag for OpenSees
    definition      % OpenSees uniaxial material definition

    K0              % Elastic stiffness
    as_Plus         % Strain hardening ratio for positive loading direction
    as_Neg          % Strain hardening ratio for negative loading direction
    Vy_Plus         % Effective yield strength for positive loading direction
    Vy_Neg          % Effective yield strength for negative loading direction (negative value)
    % Ratio of (the force at which reloading begins) to (the force corresponding to the maximum historic deformation demand) (positive loading direction)
    FprPos
    % Ratio of (the force at which reloading begins) to (the force corresponding to the absolute maximum historic deformation demand) (negative loading direction)
    FprNeg
    A_Pinch         % Ratio of reloading stiffness
    Lambda_S        % Cyclic deterioration parameter for strength deterioration
    Lambda_C        % Cyclic deterioration parameter for post-capping strength deterioration
    % Cyclic deterioration parameter for acceleration reloading stiffness deterioration
    % This is not a deterioration mode for a component with Bilinear hysteretic response.
    Lambda_A
    Lambda_K        % Cyclic deterioration parameter for unloading stiffness deterioration
    c_S             % Rate of strength deterioration. The default value is 1.0.
    c_C             % Rate of post-capping strength deterioration. The default value is 1.0.
    c_A             % Rate of accelerated reloading deterioration. The default value is 1.0.
    c_K             % Rate of unloading stiffness deterioration. The default value is 1.0.
    delta_p_Plus    % Pre-capping deflection for positive loading direction)
    delta_p_Neg     % Pre-capping deflection for negative loading direction (positive value)
    delta_pc_Plus   % Post-capping deflection for positive loading direction
    delta_pc_Neg    % Post-capping deflection for negative loading direction (positive value)
    Res_Pos         % Residual strength ratio for positive loading direction
    Res_Neg         % Residual strength ratio for negative loading direction (positive value)
    delta_u_Plus    % Ultimate deflection capacity for positive loading direction
    delta_u_Neg     % Ultimate deflection capacity for negative loading direction (positive value)
    % Rate of cyclic deterioration in the positive loading direction.
    % This parameter is used to create asymmetric hysteretic behavior for the
    % case of a composite beam. For symmetric hysteretic response use 1.0.
    D_Plus = 1.0;
    % Rate of cyclic deterioration in the negative loading direction.
    % This parameter is used to create asymmetric hysteretic behavior for the
    % case of a composite beam. For symmetric hysteretic response use 1.0.
    D_Neg = 1.0;

    includePDelta = false;              % Toggle whether or not to include implicit P-Delta effects
    stabilityCoefficient                % Stability coefficient used to calculate implicit P-Delta
    stiffnessSafety = 1.0;              % Factor on spring stiffness
    strengthSafety  = 1.0;              % Factor on spring strength
    enforceMinimumStiffness = false;    % Toggle whether to enforce a minimum stiffness between floors
    enforceMinimumStrength = false;     % Toggle whether to enforce a minimum strength between floors
    minimumRatio = 0.7;                 % Ratio used by enforceMinimumStiffness and enforceMinimumStrength
end

%=================================== Methods ==================================%
methods

% Constructor
function obj = DeterioratingSpring()

end

%---------------------------- Material definitions ----------------------------%
function matDef = bilinearMaterialDefinition(obj)
    matDef = sprintf('uniaxialMaterial Bilin %i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',...
        obj.matTag,obj.K0,obj.as_Plus,obj.as_Neg,obj.Vy_Plus,obj.Vy_Neg,obj.Lambda_S,obj.Lambda_C,obj.Lambda_A,obj.Lambda_K,...
        obj.c_S,obj.c_C,obj.c_A,obj.c_K,obj.delta_p_Plus,obj.delta_p_Neg,obj.delta_pc_Plus,obj.delta_pc_Neg,...
        obj.Res_Pos,obj.Res_Neg,obj.delta_u_Plus,obj.delta_u_Neg,obj.D_Plus,obj.D_Neg,obj.nFactor);
end

function matDef = peakOrientedMaterialDefinition(obj)
    matDef = sprintf('uniaxialMaterial ModIMKPeakOriented %i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',...
        obj.matTag,obj.K0,obj.as_Plus,obj.as_Neg,obj.Vy_Plus,obj.Vy_Neg,obj.Lambda_S,obj.Lambda_C,obj.Lambda_A,obj.Lambda_K,...
        obj.c_S,obj.c_C,obj.c_A,obj.c_K,obj.delta_p_Plus,obj.delta_p_Neg,obj.delta_pc_Plus,obj.delta_pc_Neg,...
        obj.Res_Pos,obj.Res_Neg,obj.delta_u_Plus,obj.delta_u_Neg,obj.D_Plus,obj.D_Neg);
end

function matDef = pinchingMaterialDefinition(obj)
    matDef = sprintf('uniaxialMaterial ModIMKPinching %i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',...
        obj.matTag,obj.K0,obj.as_Plus,obj.as_Neg,obj.Vy_Plus,obj.Vy_Neg,obj.FprPos,obj.FprNeg,obj.A_Pinch,obj.Lambda_S,obj.Lambda_C,obj.Lambda_A,obj.Lambda_K,...
        obj.c_S,obj.c_C,obj.c_A,obj.c_K,obj.delta_p_Plus,obj.delta_p_Neg,obj.delta_pc_Plus,obj.delta_pc_Neg,...
        obj.Res_Pos,obj.Res_Neg,obj.delta_u_Plus,obj.delta_u_Neg,obj.D_Plus,obj.D_Neg);
end

end %methods

end %classdef
