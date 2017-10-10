% test_newStrongbackDesigner.m

clear all; close all; home; %#ok<CLALL>

%================================= Definition =================================%
nStories = 6;
bldg = Strongback(nStories);

bldg.seismicDesignCategory = 'Dmax';
bldg.respModCoeff = 8;
bldg.deflAmplFact = 8;
bldg.overstrengthFactor = 3;
bldg.impFactor = 1;

%----------------------------- Units and constants ----------------------------%
bldg.units.force  = 'kip';
bldg.units.mass   = 'kslug';
bldg.units.length = 'ft';
bldg.units.time   = 'sec';

bldg.g = 32.2;       % Acceleration due to gravity

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

% --- Springs ---
spring.K0       = 1807.3;
spring.as       =  0.03;  % strain hardening ratio
spring.V_y      = 343.38;
spring.Lambda_S =  0.00;  % Cyclic deterioration parameter - strength
spring.Lambda_K =  0.00;  % Cyclic deterioration parameter - stiffness
spring.c_S      =  1.00;  % rate of deterioration - strength
spring.c_K      =  1.00;  % rate of deterioration - stiffness
spring.defl_p   =  9e99;
spring.defl_pc  =  9e99;
spring.defl_u   =  9e99;
spring.Res      =  0.30;  % residual strength ratio (relative to yield)
spring.D        =  1.00;  % rate of cyclic deterioration
spring.nFactor  =  0.00;  % elastic stiffness amplification factor

spring.definition = bilinearMaterialDefinition(1, spring);
storySpringDefinition = cell(nStories, 1);
for i = 1:nStories
    storySpringDefinition{i} = bilinearMaterialDefinition(i, spring);
end
bldg.storySpringDefinition = storySpringDefinition;

% --- Trusses ---
targetTrussDeformation = 0.002;  % Ratio of story height
bldg.storyTrussDefinition = cell(nStories,1);
trussModulus = (cumsum(bldg.storyMass,'reverse')*bldg.g)/targetTrussDeformation;
for i = 1:nStories
    bldg.storyTrussDefinition{i} = sprintf('uniaxialMaterial Elastic %i %g',i+10,trussModulus(i));
end

% --- Strongback ---
EI = 1e8;
bldg.strongbackDefinition.Area = 1;
bldg.strongbackDefinition.Modulus = sqrt(EI);
bldg.strongbackDefinition.Inertia = sqrt(EI);

%----------------------------------- Options ----------------------------------%
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;
bldg.verbose = true;

paths = fieldnames(pathOf);
for i = 1:length(paths)
    bldg.pathOf.(paths{i}) = pathOf.(paths{i});
end
bldg.pathOf.tclfunctions = '/home/petertalley/Github/OpenSees-ShearBuildings/lib';

bldg.optionsPushover.maxDrift = sum(bldg.storyHeight);
bldg.optionsPushover.test.print = 0;

bldg.optionsIDA.nMotions = 2;

gm_mat = '../ground_motions.mat';
load(gm_mat);

%================================== Run stuff =================================%

ELF = bldg.equivalentLateralForceAnalysis();

%--------------------- Plot hysteretic behavior of spring ---------------------%
anaobj = UniaxialMaterialAnalysis(spring.definition);

peakPoints  = [0 11 -11 12 -12 13 -13 14 -14 15 -15]/5;
rateType    = 'StrainRate';
rateValue   = peakPoints(2)/10;

hysteretic_pos_env = anaobj.runAnalysis([0 1.2*max(peakPoints)],rateType,rateValue);
hysteretic_neg_env = anaobj.runAnalysis([0 1.2*min(peakPoints)],rateType,rateValue);
hysteretic         = anaobj.runAnalysis(           peakPoints,  rateType,rateValue);

figure
hold on
plot(hysteretic_pos_env.disp, hysteretic_pos_env.force, 'k--')
plot(hysteretic_neg_env.disp, hysteretic_neg_env.force, 'k--')
plot(hysteretic.disp, hysteretic.force, 'r')
grid on
title('Hysteretic behavior of spring')
xlabel('Deflection (ft)')
ylabel('Force (kip)')

%---------------------------------- Pushover ----------------------------------%
F = bldg.pushoverForceDistribution();
results = bldg.pushover(F, 'TargetPostPeakRatio', 0.79);
pushover = bldg.processPushover(results, ELF);
clear results

%------------------------------ Response History ------------------------------%
IDA = bldg.incrementalDynamicAnalysis(gm_mat, pushover.periodBasedDuctility);

bldg.plotPushoverCurve(pushover)
bldg.plotPushoverDrifts(pushover)


bldg.plotIDAcurve(IDA, 'single');
xlim(gca, [0 50])
