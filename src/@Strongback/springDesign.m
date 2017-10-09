function spring = springDesign(obj,analysisResults,springGivens)
%% SPRINGDESIGN  Story spring design

spring = struct;

designStiffness = springGivens.stiffnessSafety*obj.deflAmplFact*analysisResults.storyShear./(obj.impFactor*analysisResults.allowableDrift);

designStrength  = springGivens.strengthSafety*obj.overstrengthFactor*analysisResults.storyShear;

if springGivens.enforceMinimumStiffness
    for i = 2:length(designStiffness)
        if designStiffness(i) < springGivens.minimumRatio*designStiffness(i-1)
            designStiffness(i) = springGivens.minimumRatio*designStiffness(i-1);
        end
    end
end
if springGivens.enforceMinimumStrength
    for i = 2:length(designStrength)
        if designStrength(i) < springGivens.minimumRatio*designStrength(i-1)
            designStrength(i) = springGivens.minimumRatio*designStrength(i-1);
        end
    end
end

for i = 1:obj.nStories
    if springGivens.includePDelta
        Px = sum(obj.storyMass(i:end))*obj.g;
        theta = (Px*analysisResults.allowableDrift(i)*obj.impFactor)/(analysisResults.storyShear(i)*obj.storyHeight(i)*obj.deflAmplFact);
        theta_max = min(0.5/obj.deflAmplFact,0.25);
        theta = min(theta,theta_max);
    else
        theta = 0;
    end
    spring(i).K0       = (1-theta)*designStiffness(i);  % elastic stiffness
    spring(i).as       = springGivens.as - theta;       % strain hardening ratio
    spring(i).ad       = springGivens.ad + theta;       % strain hardening ratio
    spring(i).Lambda_S = springGivens.Lambda_S;         % Cyclic deterioration parameter - strength
    spring(i).Lambda_K = springGivens.Lambda_K;         % Cyclic deterioration parameter - stiffness
    spring(i).c_S      = springGivens.c_S;              % rate of deterioration - strength
    spring(i).c_K      = springGivens.c_K;              % rate of deterioration - stiffness
    spring(i).Res      = springGivens.Res;              % residual strength ratio
    spring(i).D        = springGivens.D;                % rate of cyclic deterioration
    spring(i).nFactor  = springGivens.nFactor;          % elastic stiffness amplification factor
    spring(i).theta    = theta;

    spring(i).V_c = (1-theta)*designStrength(i);        % strength at capping
    spring(i).V_y = springGivens.C_yc*spring(i).V_c;    % effective yield strength

    spring(i).defl_y  = spring(i).V_y./spring(i).K0;                                                    % deflection at yield
    spring(i).defl_p  = (spring(i).V_c-spring(i).V_y)./((spring(i).as+theta)*spring(i).K0);             % pre-capping deflection
    spring(i).defl_pc = spring(i).V_c/(spring(i).ad*spring(i).K0);                                      % post-capping deflection
    spring(i).defl_u  = springGivens.C_upc*(spring(i).defl_y + spring(i).defl_p + spring(i).defl_pc);   % ultimate deflection capacity

    spring(i).definition = bilinearMaterialDefinition(i,spring(i));
end
end %function:springDesign
