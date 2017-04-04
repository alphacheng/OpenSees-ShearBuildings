
fpcc = 5000;
epcc = fpcc^(1/4)/4000;
Ec   = 185000*fpcc^(3/8);
rc   = fpcc/750-1.9;
xcrn = 1;
ft   = 7.5*sqrt(fpcc);
et   = ft/Ec;
rt   = 1;
xcrp = 10000;

materialDefinition = sprintf('uniaxialMaterial ConcreteCM 1 -%g -%g %g %g %g %g %g %g %g',...
    fpcc,epcc,Ec,rc,xcrn,ft,et,rt,xcrp);

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
axis manual
plot([0 0],[-10000 10000],'Color',[0.2,0.2,0.2],'LineWidth',0.25)
plot([-10000 10000],[0 0],'Color',[0.2,0.2,0.2],'LineWidth',0.25)

grid on
grid minor
