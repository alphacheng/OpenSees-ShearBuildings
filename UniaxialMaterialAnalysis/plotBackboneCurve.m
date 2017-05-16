function plotBackboneCurve(materialDefinition,endpoint,newfig)

if nargin == 2
    newfig = true;
end

anaobj = UniaxialMaterialAnalysis(materialDefinition);

rateType    = 'StrainRate';
rateValue   = 0.001;
results_pos_env = anaobj.runAnalysis([0 endpoint],rateType,rateValue);

if newfig
    figure
end
plot(results_pos_env.disp,results_pos_env.force)

end
