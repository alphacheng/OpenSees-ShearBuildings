% test_newStrongbackDesigner.m

clear all; home; %#ok<CLALL>

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

%------------------------------ Design parameters -----------------------------%

% --- Springs ---
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

% --- Trusses ---
targetTrussDeformation = 0.002;  % Ratio of story height
bldg.storyTrussDefinition = cell(nStories,1);
trussModulus = (cumsum(bldg.storyMass,'reverse')*bldg.g)/targetTrussDeformation;
for i = 1:nStories
    bldg.storyTrussDefinition{i} = sprintf('uniaxialMaterial Elastic %i %g',i+10,trussModulus(i));
end

% --- Strongback ---

%----------------------------------- Options ----------------------------------%
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;

paths = fieldnames(pathOf);
for i = 1:length(paths)
    bldg.pathOf.(paths{i}) = pathOf.(paths{i});
end
bldg.pathOf.tclfunctions = '/home/petertalley/Github/OpenSees-ShearBuildings/lib';

bldg.optionsPushover.maxDrift = sum(bldg.storyHeight);
bldg.optionsPushover.test.print = 0;

bldg.optionsIDA.nMotions = 2;

gm_mat = '../ground_motions.mat';

DEV_tol = 0.15;


%================================== Run stuff =================================%

ELF = bldg.equivalentLateralForceAnalysis();
spring = bldg.newSpringDesign(ELF, springGivens);
bldg.storySpringDefinition = spring.definition;


bldg.strongbackDefinition.Area = 1;

% Solve for appropriate EI
EI = fzero(@(x) f_dev(x, bldg, ELF, DEV_tol), 1e8);

% Do a final pass to get variables in the main namespace
bldg.strongbackDefinition.Modulus = sqrt(EI);
bldg.strongbackDefinition.Inertia = sqrt(EI);

F = bldg.pushoverForceDistribution();
results = bldg.pushover(F, 'TargetPostPeakRatio', 0.79);
pushover = bldg.processPushover(results, ELF);
VmaxI = find(pushover.baseShear == pushover.peakShear, 1);

shape = pushover.displacement_x(VmaxI, :)/pushover.displacement_x(VmaxI, end);
desired = cumsum(bldg.storyHeight)/sum(bldg.storyHeight);
deviate = max(abs( (shape - desired)./desired ));
deviation = deviate - DEV_tol;

% Plots
bldg.plotPushoverCurve(pushover)
bldg.plotPushoverDrifts(pushover)


function deviation = f_dev(EI, bldg, ELF, DEV_tol)
    bldg.strongbackDefinition.Modulus = sqrt(EI);
    bldg.strongbackDefinition.Inertia = sqrt(EI);

    F = bldg.pushoverForceDistribution();
    results = bldg.pushover(F, 'TargetPostPeakRatio', 0.95);
    pushover = bldg.processPushover(results, ELF);
    VmaxI = find(pushover.baseShear == pushover.peakShear, 1);

    shape = pushover.displacement_x(VmaxI, :)/pushover.displacement_x(VmaxI, end);
    desired = cumsum(bldg.storyHeight)/sum(bldg.storyHeight);
    deviate = max(abs( (shape - desired)./desired ));
    deviation = deviate - DEV_tol;
end
