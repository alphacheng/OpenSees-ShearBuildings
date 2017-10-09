%############################### test_R_vs_EI.m ###############################%
%                                                                              %
% Script for testing the functioning of the Strongback model.                  %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%                                                                              %
%##############################################################################%

clear all; close all; clc; %#ok<CLALL>

try

email_notify = true;
addpath('../../lib');

%################################# Definition #################################%
nStories = [3 5 7 9];
nArchetypes = length(nStories);
bldg = cellfun(@Strongback, cell(1, nArchetypes), 'UniformOutput', false);
bldg = [bldg{:}];
a = num2cell(nStories);
[bldg.nStories] = deal(a{:});
[bldg.seismicDesignCategory] = deal('Dmax');

%----------------------------- Units and constants ----------------------------%
[bldg.units] = deal(struct('force','kip','mass','kslug','length','ft','time','sec'));
% bldg.units.force = 'kip';
% bldg.units.mass  = 'kslug';
% bldg.units.length= 'ft';
% bldg.units.time  = 'sec';

[bldg.g] = deal(32.2);       % Acceleration due to gravity

[bldg.respModCoeff]       = deal(8);
[bldg.deflAmplFact]       = deal(8);
[bldg.overstrengthFactor] = deal(3);
[bldg.impFactor]          = deal(1);

%------------------------------ Story definition ------------------------------%
storyHeight = 15;
firstHeight = 20;
storyDL     = 0.080;
roofDL      = 0.030;
storyArea   = 90*90;

h = arrayfun(@(x)ones(1, x)*storyHeight, nStories, 'UniformOutput', false);
s = struct('type', '()', 'subs', {{1}});
h = cellfun(@(x)subsasgn(x, s, firstHeight), h, 'UniformOutput', false); %#ok -- works fine, not sure why CA is whining about this

[bldg.storyHeight] = deal(h{:});

d = arrayfun(@(x)ones(1, x)*storyDL, nStories, 'UniformOutput', false);
for i = 1:nArchetypes
    d{i}(end) = roofDL;
end
d = cellfun(@(x)times(x, storyArea), d, 'UniformOutput', false);
m = cellfun(@(x)rdivide(x, bldg(1).g), d, 'UniformOutput', false);

[bldg.storyMass] = deal(m{:});

%------------------------------ Spring parameters -----------------------------%
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

%------------------------------ Truss definition ------------------------------%
targetTrussDeformation = 0.002;  % Ratio of story height
storyTrussDefinition = arrayfun(@(x)cell(x, 1), nStories, 'UniformOutput', false);
trussModulus = cell(1, nArchetypes);
for i = 1:nArchetypes
    trussModulus{i} = (cumsum(bldg(i).storyMass,'reverse')*bldg(i).g)/targetTrussDeformation;
    for j = 1:nStories(i)
        storyTrussDefinition{i}{j} = sprintf('uniaxialMaterial Elastic %i %g', j+10, trussModulus{i}(j));
    end
end
[bldg.storyTrussDefinition] = deal(storyTrussDefinition{:});

strongbackDefinition.Area    = 1;
strongbackDefinition.Modulus = 1e3;
strongbackDefinition.Inertia = 1e3;

[bldg.strongbackDefinition] = deal(strongbackDefinition);

%----------------------------------- Options ----------------------------------%
[bldg.echoOpenSeesOutput] = deal(false);
[bldg.deleteFilesAfterAnalysis] = deal(true);

paths = fieldnames(pathOf);
for i = 1:nArchetypes
    for j = 1:length(paths)
        bldg(i).pathOf.(paths{j}) = pathOf.(paths{j});
    end
    bldg(i).pathOf.tclfunctions = '/home/petertalley/Github/OpenSees-ShearBuildings/lib';

    bldg(i).optionsPushover.maxDrift = sum(bldg(i).storyHeight);
    bldg(i).optionsPushover.test.print = 0;

    bldg(i).optionsIDA.nMotions = 44;
end

gm_mat = '../ground_motions.mat';

%################################# Run tests! #################################%
R_max = 8;
R_min = 1;
R_tolerance = 0.25;

% EI = logspace(0, 5.5, 1);
% EI = linspace(10^4, 10^6, 16);
EI = [142000 208000 274000 340000 406000 472000 538000 604000 670000 736000 802000 868000 934000 1000000];
R = ones(length(EI),1)*R_max;
maxMoment = zeros(length(EI), nArchetypes);
nPoints = length(EI);
for i = 1:length(EI)
    fprintf('Beginning evaluation of point %i, EI = %g %s*%s^2\n', i, EI(i), bldg(1).units.force, bldg(1).units.length)
    for j = 1:nArchetypes
        bldg(j).strongbackDefinition.Modulus = sqrt(EI(i));
        bldg(j).strongbackDefinition.Inertia = sqrt(EI(i));
    end
    R_failure = R_max;
    R_success = R_min;

    complete = false;
    failure = false;
    maxed_out = false;

    while ~complete
        fprintf('Evaluating R = %g\n', bldg(1).respModCoeff)
        ACMR   = zeros(1, nArchetypes);
        beta_t = zeros(1, nArchetypes);
        for archIndex = 1:nArchetypes
            ELF = equivalentLateralForceAnalysis(bldg(archIndex));
            spring = bldg(archIndex).springDesign(ELF,springGivens);
            bldg(archIndex).storySpringDefinition = {spring.definition}';

            F = bldg(archIndex).pushoverForceDistribution();
            pushover = bldg(archIndex).pushover(F,'TargetPostPeakRatio',0.79);
            pushover = bldg(archIndex).processPushover(pushover,ELF);

            IDA = incrementalDynamicAnalysis(bldg(archIndex), gm_mat, pushover.periodBasedDuctility);
            ACMR(archIndex)   = IDA.ACMR;
            beta_t(archIndex) = IDA.beta_t;
            r = [IDA.gm.rh];
            m = cellfun(@abs, {r.strongbackMoment}, 'UniformOutput', false);
            maxMoment(i,archIndex) = max(cell2mat(cellfun(@max, cellfun(@max, m, 'UniformOutput', false), 'UniformOutput', false)));
        end
        ACMRavg  = mean(ACMR);
        beta_avg = mean(beta_t);
        ACMR10 = FEMAP695.ACMRxx(beta_avg, 0.1);
        ACMR20 = FEMAP695.ACMRxx(beta_avg, 0.2);
        fprintf('ACMR10 = %g; ACMR20 = %g\n', ACMR10, ACMR20)
        pf = tern(ACMRavg >= ACMR10, 'passed', 'failed');
        fprintf('ACMRavg = %0g, %s (%g%% of ACMR10)\n', ACMRavg, pf, ACMRavg/ACMR10*100)
        for j = 1:nArchetypes
            pf = tern(ACMR(j) >= ACMR20, 'passed', 'failed');
            fprintf('ACMR(%i) = %0g, %s (%g%% of ACMR20)\n', j, ACMR(j), pf, ACMR(j)/ACMR20*100)
        end
        fprintf('\n')

        success = (ACMRavg >= ACMR10) && all(ACMR >= ACMR20);

        current_R = bldg(1).respModCoeff;
        if success
            if current_R == R_max
                complete = true;
                maxed_out = true;
            elseif (current_R - R_success) < R_tolerance
                complete = true;
            else
                R_success = current_R;
                [bldg.respModCoeff] = deal(R_success + (R_failure - R_success)/2);
            end
        else
            if (current_R - R_min) < R_tolerance
                complete = true;
                failure = true;
            else
                R_failure = current_R;
                [bldg.respModCoeff] = deal(R_success + (R_failure - R_success)/2);
            end
        end
    end

    if failure
        fprintf('\nAnalysis failed; R too small: %g\n', bldg(1).respModCoeff)
    else
        fprintf('\nAnalysis complete; R = %g\n', bldg(1).respModCoeff)
        fprintf('--------------------------------------------------------------------------------\n')
    end

    if maxed_out
        fprintf('Maxed out R; ACMR only goes up from here\n')
        break
    end

    R(i) = bldg(1).respModCoeff;
end

% --Resuming from previous operation--
load prev.mat
EI = linspace(10^4, 10^6, 16);
R = vertcat(R_save, R);
maxMoment = vertcat(mM_save, maxMoment);

% Plot R vs EI
plot_handles(1) = figure;
semilogx(EI, R, '*-')
ylim([0 10])
ylabel('Response Modification Coefficient')
xlabel('EI of Strongback (kip*ft^2)')
grid on
title(sprintf('R vs EI -- %i building test suite, %i ground motions', nArchetypes, bldg(1).optionsIDA.nMotions))

% Plot max moment vs EI
plot_handles(2) = figure;
plot(EI, maxMoment, '*')
ylabel('Absolute maximum moment in strongback')
xlabel('EI of Strongback (kip*ft^2)')
grid on

% Plot max shear vs EI


% Email me the plots when it's done
if email_notify
    plot_filefolder = pwd;
    plot_filename   = sprintf('R_vs_EI_%g.fig', now);
    plot_filepath   = fullfile(plot_filefolder, plot_filename);
    savefig(plot_handles, plot_filepath, 'compact');
    if filesize(plot_filepath) > Email.max_attachment_size('gmail')
        % if the plot file is bigger than what gmail allows, need to send a dropbox link instead
        dropfile = fullfile(pathOf.Dropbox, 'figure_share', [plot_filename '.fig']);
        movefile(plot_filepath, dropfile);
        link = Dropbox.sharelink(dropfile);
        Email.mutt_send('R vs EI plot finished', Email.contacts('utk'), '', link);
    else
        Email.mutt_send('R vs EI plot finished', Email.contacts('utk'), plot_filepath);
    end
end

catch err
    save('err.mat', 'err')
    Email.mutt_send('An error occurred', Email.contacts('utk'), 'err.mat', err.message);
end
