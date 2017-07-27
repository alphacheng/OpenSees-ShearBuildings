% Test the workings of the mdofShearTrussBuilding class.
clear all;close all;clc; %#ok<CLALL>

%##############################################################################%
%% Define Building
nStories = 3;
bldg = mdofShearTrussBuilding(nStories);

% Units
bldg.units.force = 'kip';
bldg.units.mass  = 'kslug';
bldg.units.length= 'ft';
bldg.units.time  = 'sec';

bldg.g = 32.2;       % Acceleration due to gravity (ft/s^2)

storyHeight = 15;    % Story height (ft)
firstHeight = 20;    %
storyDL     = 0.080; % Story dead loads (ksf)
roofDL      = 0.030; % Roof dead load (ksf)
storyArea   = 90*90; % Story areas (ft^2)

bldg.storyHeight    = ones(1,nStories)*storyHeight;
bldg.storyHeight(1) = firstHeight;

storyDL        = ones(1,nStories)*storyDL;
storyDL(end)   = roofDL;
bldg.storyMass = (storyDL*storyArea)/bldg.g;

bldg.seismicDesignCategory = 'Dmax';
bldg.respModCoeff = 8;
bldg.deflAmplFact = 8;
bldg.overstrengthFactor = 3;
bldg.impFactor = 1;

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

%##############################################################################%
%% Analysis Options
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = true;

bldg.verbose   = true ; % Toggle verbose output
paths = fieldnames(pathOf);
for i = 1:length(paths)
    bldg.pathOf.(paths{i}) = pathOf.(paths{i});
end
bldg.pathOf.tclfunctions = '/home/petertalley/Github/OpenSees-ShearBuildings/lib';

%-------------------------------------------------------------------------------
% Pushover analysis options
bldg.optionsPushover.controlStory = 'roof';
bldg.optionsPushover.stepSize = 0.001;
bldg.optionsPushover.maxDrift = min(bldg.storyHeight);

bldg.optionsPushover.constraints.type = 'Plain';
bldg.optionsPushover.constraints.penalty.alphaS = 1.0e12;
bldg.optionsPushover.constraints.penalty.alphaM = 1.0e12;

bldg.optionsPushover.test.type       = 'NormDispIncr';
bldg.optionsPushover.test.tolerance  = [1e-5,1e-4,1e-3];
bldg.optionsPushover.test.iterations = 10;
bldg.optionsPushover.test.print      = 0;
bldg.optionsPushover.test.normType   = 2;

bldg.optionsPushover.algorithm = { 'KrylovNewton','Newton','ModifiedNewton' };

%-------------------------------------------------------------------------------
% Response history options
bldg.damping_ModeA  = 1;             % Mode A for rayleigh damping
bldg.damping_ModeB  = 3;             % Mode B for rayleigh damping
bldg.damping_RatioA = 0.05;          % Damping ratio for mode A
bldg.damping_RatioB = 0.05;          % Damping ratio for mode B

bldg.optionsResponseHistory.constraints.type = 'Transformation';
bldg.optionsResponseHistory.constraints.penalty.alphaS = 1.0e12;
bldg.optionsResponseHistory.constraints.penalty.alphaM = 1.0e12;

bldg.optionsResponseHistory.test.type       = 'NormDispIncr';
bldg.optionsResponseHistory.test.tolerance  = [1e-5,1e-4,1e-3];
bldg.optionsResponseHistory.test.iterations = 10;
bldg.optionsResponseHistory.test.print      = 0;
bldg.optionsResponseHistory.test.normType   = 2;

bldg.optionsResponseHistory.algorithm = { 'KrylovNewton','Newton','ModifiedNewton' };

%-------------------------------------------------------------------------------
% Incremental dynamic analysis options
bldg.optionsIDA.nMotions = 44;                              % Number of ground motions to analyze
bldg.optionsIDA.tExtra = 5;                                 % Extra analysis time after end of ground motion
bldg.optionsIDA.collapseDriftRatio = 0.05;                  % Story drift ratio that defines collapse
bldg.optionsIDA.ST = [0.5:0.5:2.5 2.75:0.25:8];
bldg.optionsIDA.rating_DR  = 'C';
bldg.optionsIDA.rating_TD  = 'C';
bldg.optionsIDA.rating_MDL = 'C';
bldg.optionsIDA.shortCircuit = true;
bldg.optionsIDA.ST_tol = 0.01;
bldg.optionsIDA.ST_step = 0.5;

%##############################################################################%
%% Run Test

nHistories = length(bldg.optionsIDA.ST);
gm_mat = '../ground_motions.mat';

ELF = bldg.ELFanalysis();
spring = bldg.springDesign(ELF,springGivens);
bldg.storySpringDefinition = {spring.definition}';
bldg.storyTrussDefinition = cell(nStories,1);
targetTrussDeformation = 0.002;  % Ratio of story height
trussModulus = (cumsum(bldg.storyMass,'reverse')*bldg.g)/targetTrussDeformation;
for i = 1:nStories
    bldg.storyTrussDefinition{i} = sprintf('uniaxialMaterial Elastic %i %g',i*10,trussModulus(i));
end

% Standard setup - P-Delta only modelled explicitly
F = bldg.pushoverForceDistribution();
push = bldg.pushover(F,'TargetPostPeakRatio',0.75);
push = bldg.processPushover(push,ELF);
bldg.plotPushoverCurve(push);
bldg.plotPushoverDrifts(push)

IDA = bldg.incrementalDynamicAnalysis(gm_mat,push.periodBasedDuctility);
for plotMode = {'single','multiple'}
    plotIDAcurve(bldg,IDA,plotMode{1})
end
