%################################ PDeltaTest.m ################################%
%                                                                              %
% Script for testing the weirdness involved with explicit P-Delta effects.     %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%##############################################################################%

clear all; close all; clc; %#ok<CLALL>

%################################# Definition #################################%
nStories = 3;
bldg = mdofShearTrussBuilding(nStories);
bldg.seismicDesignCategory = 'Dmax';
bldg.respModCoeff = 8;
bldg.deflAmplFact = 8;
bldg.overstrengthFactor = 3;
bldg.impFactor = 1;

%----------------------------- Units and constants ----------------------------%
bldg.units.force = 'kip';
bldg.units.mass  = 'kslug';
bldg.units.length= 'ft';
bldg.units.time  = 'sec';

bldg.g = 32.2;       % Acceleration due to gravity (ft/s^2)

%------------------------------ Story definition ------------------------------%
storyHeight = 15;
firstHeight = 20;
storyDL     = 0.080;
roofDL      = 0.030;
storyArea   = 90*90;

bldg.storyHeight    = ones(1,nStories)*storyHeight;
bldg.storyHeight(1) = firstHeight;

storyDL        = ones(1,nStories)*storyDL;
storyDL(end)   = roofDL;
bldg.storyMass = (storyDL*storyArea)/bldg.g;

%----------------------------- Spring information -----------------------------%
springGivens.as       =  0.03;  % strain hardening ratio
springGivens.Lambda_S = 10.00;  % Cyclic deterioration parameter - strength
springGivens.Lambda_K = 10.00;  % Cyclic deterioration parameter - stiffness
springGivens.c_S      =  1.00;  % rate of deterioration - strength
springGivens.c_K      =  1.00;  % rate of deterioration - stiffness
springGivens.Res      =  0.30;  % residual strength ratio (relative to yield)
springGivens.D        =  1.00;  % rate of cyclic deterioration
springGivens.nFactor  =  0.00;  % elastic stiffness amplification factor
springGivens.C_yc     =  0.80;  % ratio of yield strength to capping strength
springGivens.C_upc    = 20.00;  % ratio of ultimate deflection to u_y + u_p + u_pc
springGivens.ad       =  0.10;  % deterioration stiffness ratio -- higher values mean faster deterioration
springGivens.includePDelta = false;

springGivens.stiffnessSafety = 1.0;
springGivens.strengthSafety  = 1.0;

springGivens.enforceMinimumStiffness = false;
springGivens.enforceMinimumStrength = false;
springGivens.minimumRatio = 0.7;

targetTrussDeformation = 0.005;  % Ratio of story height

%------------------------------ Analysis options ------------------------------%
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;
bldg.verbose   = true ; % Toggle verbose output

paths = fieldnames(pathOf);
for i = 1:length(paths)
    bldg.pathOf.(paths{i}) = pathOf.(paths{i});
end
bldg.pathOf.tclfunctions = '/home/petertalley/Github/OpenSees-ShearBuildings/lib';

%################################### Design ###################################%
ELF = bldg.ELFanalysis();
spring = bldg.springDesign(ELF,springGivens);
bldg.storySpringDefinition = {spring.definition}';

bldg.storyTrussDefinition = cell(nStories,1);
trussModulus = (cumsum(bldg.storyMass,'reverse')*bldg.g)/targetTrussDeformation;
for i = 1:nStories
    bldg.storyTrussDefinition{i} = sprintf('uniaxialMaterial Elastic %i %g',i*10,trussModulus(i));
end

%#################################### Tests ###################################%

%------------------------------ P-Delta included ------------------------------%
F = bldg.pushoverForceDistribution();
results1 = pushover(bldg,F,'TargetPostPeakRatio',0.75);
results1 = processPushover(bldg,results1,ELF);
fig = bldg.plotPushoverCurve(results1);

%----------------------------- No P-Delta included ----------------------------%
bldg.includeExplicitPDelta = false;
F = bldg.pushoverForceDistribution();
results2 = pushover(bldg,F,'TargetPostPeakRatio',0.75);
results2 = processPushover(bldg,results2,ELF);
fig = bldg.plotPushoverCurve(results2,fig);

legend('With P-Delta','No P-Delta','Location','southeast')
