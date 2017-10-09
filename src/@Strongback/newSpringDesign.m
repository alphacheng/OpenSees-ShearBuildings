function spring = newSpringDesign(obj,ELF,springGivens)
%% SPRINGDESIGN  Story spring design

% Alias things to shorter titles
F_k = springGivens.stiffnessSafety;
F_f = springGivens.strengthSafety;

I = obj.impFactor;
C_d = obj.deflAmplFact;

avgShear = ELF.baseShear/obj.nStories;
avgDrift = mean(ELF.allowableDrift);

% Initialize results
spring = struct;

designStiffness = F_k*C_d*avgShear./(I*avgDrift);
designStrength  = F_f*obj.overstrengthFactor*avgShear;

spring.K0       = designStiffness;  % elastic stiffness
spring.as       = springGivens.as;       % strain hardening ratio
spring.ad       = springGivens.ad;       % strain hardening ratio
spring.Lambda_S = springGivens.Lambda_S;         % Cyclic deterioration parameter - strength
spring.Lambda_K = springGivens.Lambda_K;         % Cyclic deterioration parameter - stiffness
spring.c_S      = springGivens.c_S;              % rate of deterioration - strength
spring.c_K      = springGivens.c_K;              % rate of deterioration - stiffness
spring.Res      = springGivens.Res;              % residual strength ratio
spring.D        = springGivens.D;                % rate of cyclic deterioration
spring.nFactor  = springGivens.nFactor;          % elastic stiffness amplification factor

spring.V_c = designStrength;        % strength at capping
spring.V_y = springGivens.C_yc*spring.V_c;    % effective yield strength

spring.defl_y  = spring.V_y./spring.K0;                                                    % deflection at yield
spring.defl_p  = (spring.V_c-spring.V_y)./(spring.as*spring.K0);             % pre-capping deflection
spring.defl_pc = spring.V_c/(spring.ad*spring.K0);                                      % post-capping deflection
spring.defl_u  = springGivens.C_upc*(spring.defl_y + spring.defl_p + spring.defl_pc);   % ultimate deflection capacity

definition = cell(obj.nStories, 1);
for i = 1:obj.nStories
    definition{i} = bilinearMaterialDefinition(i,spring);
end
spring.definition = definition;

end %function:springDesign
