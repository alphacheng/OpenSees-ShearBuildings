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

%------------------------------ Design parameters -----------------------------%

% --- Springs ---
springGivens.as       =  0.03;  % strain hardening ratio
springGivens.Lambda_S =  0.00;  % Cyclic deterioration parameter - strength
springGivens.Lambda_K =  0.00;  % Cyclic deterioration parameter - stiffness
springGivens.c_S      =  1.00;  % rate of deterioration - strength
springGivens.c_K      =  1.00;  % rate of deterioration - stiffness
springGivens.Res      =  1.25;  % residual strength ratio (relative to yield)
springGivens.D        =  1.00;  % rate of cyclic deterioration
springGivens.nFactor  =  0.00;  % elastic stiffness amplification factor
springGivens.C_yc     =  0.80;  % ratio of yield strength to capping strength
springGivens.C_upc    = 20.00;  % ratio of ultimate deflection to u_y + u_p + u_pc
springGivens.ad       =  0.10;  % deterioration stiffness ratio -- higher values mean faster deterioration

springGivens.stiffnessSafety = 1;
springGivens.strengthSafety  = 2;

ELF = bldg.equivalentLateralForceAnalysis();
spring = bldg.newSpringDesign(ELF, springGivens);
bldg.storySpringDefinition = spring.definition;

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

bldg.optionsIDA.nMotions = 1;

gm_mat = '../ground_motions.mat';
load(gm_mat);

use_existing_results = false;
existing_results = './interesting_results.mat';

%================================== Run stuff =================================%

if use_existing_results
    load(existing_results)
    load(gm_mat)
else
    %-------------------------------- Pushover --------------------------------%
    F = bldg.pushoverForceDistribution();
    results = bldg.pushover(F, 'TargetPostPeakRatio', 0.79);
    pushover = bldg.processPushover(results, ELF);
    clear results

    %---------------------------- Response History ----------------------------%
    IDA = bldg.incrementalDynamicAnalysis(gm_mat, pushover.periodBasedDuctility);
end

%------------------------------- Standard plots -------------------------------%
fig = bldg.plotPushoverCurve(pushover);
ax = fig.Children;
hold(ax, 'on')
strongbackBaseShear = sum(pushover.strongbackShear, 2);
plot(ax, pushover.roofDrift, strongbackBaseShear)
plot(ax, xlim(ax), [ELF.baseShear, ELF.baseShear], '--')
ylabel('Force (kip)')
legend(ax, 'Total base shear', 'Strongback base shear', 'ELF base shear')

%---------------- Ratio of base shear to strongback base shear ----------------%
% Pushover
figure
ratio = pushover.baseShear ./ strongbackBaseShear;
plot(pushover.roofDrift, ratio);
grid on
xlabel('Roof drift (ft)')
ylabel('Ratio of base shear to strongback base shear')

%------------------------------ Normalized shape ------------------------------%
bldg.plotPushoverDrifts(pushover)

VmaxI = find(pushover.baseShear == pushover.peakShear, 1);
shape = pushover.displacement_x(VmaxI, :)/pushover.displacement_x(VmaxI, end);
desired = cumsum(bldg.storyHeight)/sum(bldg.storyHeight);

figure
plot([0 desired], [0 cumsum(bldg.storyHeight)], '*-')
hold on
plot([0 shape], [0 cumsum(bldg.storyHeight)], '*-')
grid on
xlim([0 1.2])
ylabel('Height (ft)')
xlabel('Normalized displacement \delta/\delta_r')
xlim([0 1])
title('Deflected shapes')
legend('Desired (no damage concentration)', 'Actual')

%------------------- Strongback V & M envelopes (pushover) --------------------%
% Create moment plot figure
figure
ax_m = gca;
hold on
grid on
% Create shear plot figure
figure
ax_s = gca;
hold on
grid on

% Plot moment envelope
maxMoment = max(abs(pushover.strongbackMoment));
plot(ax_m, maxMoment, [0 cumsum(bldg.storyHeight)], '*-')

% Plot shear envelope
maxShear = max(abs(pushover.strongbackShear));
plot(ax_s, maxShear, [0 cumsum(bldg.storyHeight)], '*-')

title(ax_m, sprintf('Moment envelope, %i-story building\n EI=%g %s-%s^2, pushover',...
    bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia,...
    bldg.units.force, bldg.units.length))
title(ax_s, sprintf('Shear envelope, %i-story building\n EI=%g %s-%s^2, pushover',...
    bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia,...
    bldg.units.force, bldg.units.length))

xlabel(ax_m, sprintf('Moment (%s-%s)', bldg.units.force, bldg.units.length))
xlabel(ax_s, sprintf('Shear (%s)', bldg.units.force))
ylabel(ax_m, sprintf('Height (%s)', bldg.units.length))
ylabel(ax_s, sprintf('Height (%s)', bldg.units.length))


%------------------------------ IDA drift ratios ------------------------------%
bldg.plotIDAcurve(IDA, 'single')
title(gca, 'Incremental Dynamic Analysis')

figure
hold on
legend_entries = cell(1, bldg.optionsIDA.nMotions);
for i = 1:bldg.optionsIDA.nMotions
    r = IDA.gm(i).rh(IDA.gm(i).collapseIndex-1);
    storyDriftRatio = r.storyDrift ./ bldg.storyHeight;
    [row, ~] = find(abs(storyDriftRatio) == r.maxDriftRatio);
    plot(abs(storyDriftRatio(row, :))*100, cumsum(bldg.storyHeight))
    legend_entries{i} = IDA.gm(i).ID;
end
legend(legend_entries)
ylim([0 sum(bldg.storyHeight)+bldg.storyHeight(1)])
grid on
xlabel('Story drift ratio (%)')
ylabel('Height (ft)')
title('Story drift ratios at intensity immediately prior to collapse')


%---------------------- Strongback V & M envelopes (IDA) ----------------------%
% Create moment plot figure
figure
ax_m = gca;
hold on
grid on
% Create shear plot figure
figure
ax_s = gca;
hold on
grid on

for i = 1:bldg.optionsIDA.nMotions
    % Plot moment envelope
    maxMoment = max(abs(IDA.gm(i).rh(IDA.gm(i).collapseIndex-1).strongbackMoment));
    plot(ax_m, maxMoment, [0 cumsum(bldg.storyHeight)], '*-')

    % Plot shear envelope
    maxShear = max(abs(IDA.gm(i).rh(IDA.gm(i).collapseIndex-1).strongbackShear));
    plot(ax_s, maxShear, [0 cumsum(bldg.storyHeight)], '*-')
end
title(ax_m, sprintf('Moment envelope, %i-story building\n EI=%g %s-%s^2, gm %s-%s',...
    bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia, bldg.units.force,...
    bldg.units.length, ground_motions(1).ID, ground_motions(bldg.optionsIDA.nMotions).ID))
title(ax_s, sprintf('Shear envelope, %i-story building\n EI=%g %s-%s^2, gm %s-%s',...
    bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia, bldg.units.force,...
    bldg.units.length, ground_motions(1).ID, ground_motions(bldg.optionsIDA.nMotions).ID))
legend(ax_m, {ground_motions(1:bldg.optionsIDA.nMotions).ID})
legend(ax_s, {ground_motions(1:bldg.optionsIDA.nMotions).ID})

xlabel(ax_m, sprintf('Moment (%s-%s)', bldg.units.force, bldg.units.length))
xlabel(ax_s, sprintf('Shear (%s)', bldg.units.force))
ylabel(ax_m, sprintf('Height (%s)', bldg.units.length))
ylabel(ax_s, sprintf('Height (%s)', bldg.units.length))
