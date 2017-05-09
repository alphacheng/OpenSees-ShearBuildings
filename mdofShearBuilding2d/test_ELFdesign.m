% test_ELFdesign.m
% Units: kslug, kip, ft, sec

%##############################################################################%
%% Startup

tic
clear; close all; clc;
neededPaths={'../'; ...
             '../UniaxialMaterialAnalysis'};

for i = 1:length(neededPaths)
    addpath(neededPaths{i})
end

%##############################################################################%
%% Define Building
nStories = 3;
bldg = mdofShearBuilding2d(nStories);

bldg.g = 32.2;                                  % Acceleration due to gravity (ft/s^2)

bldg.storyHeight = [20 15 15];                  % Story heights (ft)
storyDL = [0.080 0.080 0.030];                  % Story dead loads (ksf)
storyArea = 90*90*ones(1,nStories);             % Story areas (ft^2)
bldg.storyMass = (storyDL .* storyArea)/bldg.g; % Story masses (kslug)

bldg.seismicDesignCategory = 'Dmax';
bldg.respModCoeff = 8;
bldg.deflAmplFact = 8;
bldg.overstrengthFactor = 3;
bldg.impFactor = 1;

springGivens.as      = 0.04;    % strain hardening ratio
springGivens.Lambda  = 8.00;    % Cyclic deterioration parameter
springGivens.c       = 1.00;    % rate of deterioration
springGivens.Res     = 0.30;    % residual strength ratio
springGivens.D       = 1.00;    % rate of cyclic deterioration
springGivens.nFactor = 0.00;    % elastic stiffness amplification factor
springGivens.C_yc    = 0.90;    % ratio of yield strength to capping strength
springGivens.C_pcp   = 1.00;    % ratio of post-capping deflection to pre-capping deflection
springGivens.C_upc   = 20.0;    % ratio of ultimate deflection to u_y + u_p + u_pc

% Units -- used for descriptions only
units.force = 'kip';
units.mass  = 'kslug';
units.length= 'ft';
units.time  = 'sec';

%##############################################################################%
%% Analysis Options
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = true;

verbose     = true ;    % Toggle verbose output
runPushover = true ;    % Toggle pushover analysis
runIDA      = false;    % Toggle IDA

% Equivalent lateral force options
iterate = false;             % Select whether to do iteration
iterOption = 'period';      % Variable to use for convergence: 'period' or 'overstrength'

diffTol = 1e-3;             % Tolerance for iteration

% Pushover analysis options
bldg.pushover_stepSize   = 0.001;
bldg.pushover_maxDrift   = 100;

% Incremental dynamic analysis options
nMotions = 6;                              % Number of ground motions to analyze
SF2 = [0:0.25:1.5 , 2:0.5:5 , 5.75:0.75:8]; % Scale factors to use for each IDA curve

%##############################################################################%
%% Equivalent lateral force procedure

if verbose; fprintf('Running equivalent lateral force procedure...\n'); end

iterating = true;
currDiff = 1+diffTol;
while iterating == true
    %% Define springs
    resultsELF = bldg.ELFanalysis();
    spring = bldg.springDesign(resultsELF,springGivens);

    bldg.storySpringDefinition = spring.definition;

    if ~iterate
        iterating = false;
    else
        switch lower(iterOption)
            case 'period'
                eigenvals = bldg.eigenvalues();
                calculatedPeriod = (sqrt(eigenvals(1))/(2*pi))^-1;
                prevDiff = currDiff;
                currDiff = bldg.fundamentalPeriod - calculatedPeriod;
                bldg.fundamentalPeriod = calculatedPeriod;
                if prevDiff == currDiff || abs(currDiff) < diffTol
                    iterating = false;
                end
            case 'overstrength'
                [~,eigenvecs] = bldg.eigenvalues;
                F = bldg.storyMass' .* eigenvecs;
                resultsPushover = bldg.pushover(F,'TargetPostPeakRatio',0.75);
                Vmax = max(resultsPushover.baseShear);
                calcOverstr = Vmax/resultsELF.baseShear;
                prevDiff = currDiff;
                currDiff = bldg.overstrengthFactor - calcOverstr;
                bldg.overstrengthFactor = calcOverstr;
                if prevDiff == currDiff || abs(currDiff) < diffTol
                    iterating = false;
                end
        end
    end
end

%##############################################################################%
%% Pushover analysis
if runPushover
if verbose
    fprintf('Running pushover analysis...\n');
    pushover_tic = tic;
end

[~,eigenvecs] = bldg.eigenvalues;
F = bldg.storyMass' .* eigenvecs;

resultsPushover = bldg.pushover(F,'TargetPostPeakRatio',0.75);

peakShear = max(resultsPushover.baseShear);
peakShearIndex = resultsPushover.baseShear == peakShear;

postPeakIndex = resultsPushover.roofDrift > resultsPushover.roofDrift(peakShearIndex);
postPeakShear = resultsPushover.baseShear(postPeakIndex);
postPeakDrift = resultsPushover.totalDrift(postPeakIndex,:);

peakShear80 = 0.8*peakShear;
peakShear80Drift = interp1(postPeakShear,postPeakDrift,peakShear80);

peakStoryDrift = resultsPushover.storyDrift(peakShearIndex,:);
postPeakStoryDrift = resultsPushover.storyDrift(postPeakIndex,:);
peak80StoryDrift = interp1(postPeakShear,postPeakStoryDrift,peakShear80);

peakStoryDriftRatio   = peakStoryDrift./bldg.storyHeight;
peak80StoryDriftRatio = peak80StoryDrift./bldg.storyHeight;

fprintf('80%% Peak Shear = %g %s\n',peakShear80,units.force);
fprintf('Roof drift at 80%% Peak Shear = %g %s\n',peakShear80Drift(bldg.nStories),units.length);

figure
hold on
plot(resultsPushover.roofDrift,resultsPushover.baseShear,'-')
grid on
grid minor
xlabel(sprintf('Roof drift (%s)',units.length))
ylabel(sprintf('Base shear (%s)',units.force))
title('Pushover analysis')

figure
subplot(1,2,1)
barh(F/sum(F),0.1)
grid on
title('Pushover force distribution')
ylabel('Story')
subplot(1,2,2)
hold on
plot(peakStoryDriftRatio  ,1:bldg.nStories,'*-')
plot(peak80StoryDriftRatio,1:bldg.nStories,'*-')
grid on
grid minor
ylim([0.5 bldg.nStories+0.5])
ylabel('Story')
xlabel('Story Drift Ratio')
title('Pushover story drifts')
legend('V_{max}','V_{80}','Location','Southeast')

if verbose
    pushover_time = toc(pushover_tic);
    fprintf('Pushover analysis took %g seconds.\n',pushover_time);
end

end

%##############################################################################%
%% Incremental Dynamic Analysis
if runIDA

if verbose
    fprintf('Running incremental dynamic analysis...\n');
    ida_tic = tic;
end

load('ground_motions.mat');
SMT = FEMAP695_SMT(bldg.fundamentalPeriod,bldg.seismicDesignCategory);
ST  = SMT*SF2;
figure
hold on
legendentries = cell(nMotions,1);
resultsIDA = cell(nMotions,length(SF2));
for i = 1:nMotions
    gmfile = bldg.scratchFile(sprintf('acc%s.acc',ground_motions(i).ID));
    gmfid = fopen(gmfile,'w+');
    for k = 1:ground_motions(i).numPoints
        fprintf(gmfid,'%f\n',ground_motions(i).normalized_acceleration(k)*bldg.g);
    end
    fclose(gmfid);
    dt      = ground_motions(i).dt;
    SF1     = FEMAP695_SF1(bldg.fundamentalPeriod,bldg.seismicDesignCategory);
    tend    = max(ground_motions(i).time) + 5;

    % Vary scale factor
    maxDriftRatio = zeros(length(SF2),1);
    parfor j = 1:length(SF2)
        if verbose
            fprintf('Calculating IDA point for %s, SF2 = %g\n',ground_motions(i).ID,SF2(j));
        end
        SF = SF1*SF2(j);
        resultsIDA{i,j} = bldg.responseHistory(gmfile,dt,SF,tend,ground_motions(i).ID,j);
        maxDriftRatio(j) = max(max(abs(resultsIDA{i,j}.storyDrift))./bldg.storyHeight);
    end

    plot(maxDriftRatio,ST,'o-')
    legendentries{i} = ground_motions(i).ID;

    if bldg.deleteFilesAfterAnalysis
        delete(gmfile)
    end
end

grid on
xlabel('Maximum story drift ratio')
ylabel('Ground motion intensity, S_T (g)')
legend(legendentries)

% Plot sample response history
plotSampleResponse(resultsIDA{1,5})

if verbose
    ida_time = toc(ida_tic);
    fprintf('Incremental dynamic analysis took %g seconds.\n',ida_time);
end

end

%##############################################################################%
%% Plot backbone curves
figure
hold on
backbone_ax = cell(nStories,1);
for i = 1:nStories
    backbone_ax{i} = subplot(nStories,1,i);
    materialDefinition = bldg.storySpringDefinition{i};
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';
    plotBackboneCurve(materialDefinition,spring.defl_u(i),false)
    xlim([0 1.1*(spring.defl_y(i)+spring.defl_p(i)+spring.defl_pc(i))])
    title(sprintf('Backbone curve for story %i',i))
    xlabel(sprintf('Deflection (%s)',units.length))
    ylabel(sprintf('Force (%s)',units.force))
    grid on
end

%##############################################################################%
%% Cleanup

for i = 1:length(neededPaths)
    rmpath(neededPaths{i})
end

toc
