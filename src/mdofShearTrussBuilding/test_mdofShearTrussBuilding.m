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
bldg.optionsIDA.nMotions = 12;                              % Number of ground motions to analyze
bldg.optionsIDA.tExtra = 5;                                 % Extra analysis time after end of ground motion
bldg.optionsIDA.collapseDriftRatio = 0.05;                  % Story drift ratio that defines collapse
bldg.optionsIDA.ST = [0.5:0.5:2.5 2.75:0.25:8];
bldg.optionsIDA.rating_DR  = 'C';
bldg.optionsIDA.rating_TD  = 'C';
bldg.optionsIDA.rating_MDL = 'C';
bldg.optionsIDA.shortCircuit = true;

%##############################################################################%
%% Run Test

nHistories = length(bldg.optionsIDA.ST);

ELF = bldg.ELFanalysis();
spring = bldg.springDesign(ELF,springGivens);
bldg.storySpringDefinition = {spring.definition}';
bldg.storyTrussDefinition = cell(nStories,1);
a_t = 100.0;
trussModulus = a_t*[spring.K0].*bldg.storyHeight;
for i = 1:nStories
    bldg.storyTrussDefinition{i} = sprintf('uniaxialMaterial Elastic %i %g',i*10,trussModulus(i));
end

SF1 = FEMAP695_SF1(bldg.fundamentalPeriod,bldg.seismicDesignCategory);
SF2 = bldg.optionsIDA.ST/FEMAP695_SMT(bldg.fundamentalPeriod,bldg.seismicDesignCategory);
SF  = SF1*SF2;

load('../ground_motions.mat')
gm(bldg.optionsIDA.nMotions) = struct;
parfor gmIndex = 1:bldg.optionsIDA.nMotions
    gmID = ground_motions(gmIndex).ID;
    gmFile = scratchFile(bldg,sprintf('acc%s.acc',gmID));
    dlmwrite(gmFile,ground_motions(gmIndex).normalized_acceleration*bldg.g);
    dt     = ground_motions(gmIndex).dt;
    tEnd   = ground_motions(gmIndex).time(end);

    rh = cell(nHistories,1); % Operating on structs in isolate simplifies adding things to them
    for index = 1:nHistories
        rh{index} = bldg.responseHistory(gmFile,dt,SF(index),tEnd,gmID,index); %#ok<PFBNS>
        rh{index}.energy = bldg.energyCriterion(rh{index});
        rh{index}.maxDriftRatio = max(max(abs(rh{index}.storyDrift))./bldg.storyHeight);
    end
    rh = [rh{:}];
    temp = [rh.energy];
    gm(gmIndex).maxDriftRatio = [rh.maxDriftRatio];
    gm(gmIndex).collapseIndex = find(~isnan([temp.failureIndex]),1);
    gm(gmIndex).rh = rh;

    if bldg.deleteFilesAfterAnalysis
        delete(gmFile)
    end
end

figure
hold on
for gmIndex = 1:bldg.optionsIDA.nMotions
    collapse = gm(gmIndex).collapseIndex;
    plot([0 gm(gmIndex).maxDriftRatio(1:collapse)*100],[0 bldg.optionsIDA.ST(1:collapse)],'ko-')
    plot(gm(gmIndex).maxDriftRatio(collapse)*100,bldg.optionsIDA.ST(collapse),'k.','MarkerSize',10)
end
xlim([0 3*bldg.optionsIDA.collapseDriftRatio*100])
xlabel('Maximum Story Drift Ratio (%)')
ylabel('Ground Motion Intensity S_T (g)')
grid on

% for i = [1 collapseIndex-1 collapseIndex nHistories]
%     figure
%     plot(gm(gmIndex).rh(i).time,[gm(gmIndex).rh(i).energy.earthquake gm(gmIndex).rh(i).energy.norm_gravity])
%     yl = ylim;
%     ylim([0 yl(2)])
%     grid on
%     xlabel(sprintf('Time (%s)',bldg.units.time))
%     ylabel(sprintf('Energy (%s*%s)',bldg.units.length,bldg.units.force))
% end
