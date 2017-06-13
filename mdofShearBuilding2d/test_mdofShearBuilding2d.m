% test_ELFdesign.m
% Units: kslug, kip, ft, sec

%##############################################################################%
%% Startup

tic
clear all; close all; clc;   %#ok<CLALL>
neededPaths={'../'
             '../UniaxialMaterialAnalysis'};

for i = 1:length(neededPaths)
    addpath(neededPaths{i})
end

results = struct;

test_options
if isempty(gcp)
    parpool;
end

%##############################################################################%
%% Equivalent lateral force procedure

if bldg.verbose; fprintf('Running equivalent lateral force procedure...\n'); end

iterating = true;
currDiff = 1+diffTol;
nIter = 0;
while iterating == true
    nIter = nIter + 1;
    %% Define springs
    results.ELF = bldg.ELFanalysis();
    spring = bldg.springDesign(results.ELF,springGivens);

    storySpringDefinition = cell(nStories,1);
    for i = 1:nStories
        storySpringDefinition{i} = spring(i).definition;
    end
    bldg.storySpringDefinition = storySpringDefinition;

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
                    fprintf('Did %i iterations.\n',nIter)
                end
            case 'overstrength'
                [~,eigenvecs] = bldg.eigenvalues;
                F = bldg.storyMass' .* eigenvecs;
                results.pushover = bldg.pushover(F,'TargetPostPeakRatio',0.75);
                Vmax = max(results.pushover.baseShear);
                calcOverstr = Vmax/results.ELF.baseShear;
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
if bldg.verbose
    fprintf('Running pushover analysis...\n');
    pushover_tic = tic;
end

[~,eigenvecs] = bldg.eigenvalues;

F = bldg.storyMass' .* eigenvecs;
if min(F) < 0
    F = -F;
end

results.pushover = bldg.pushover(F,'TargetPostPeakRatio',0.75);

figure
hold on
plot(results.pushover.roofDrift,results.pushover.baseShear,'-')
grid on
grid minor
xlabel(sprintf('Roof drift (%s)',bldg.units.length))
ylabel(sprintf('Base shear (%s)',bldg.units.force))
title('Pushover analysis')

switch results.pushover.exitStatus
case 'Analysis Successful'
    calcOverstr = results.pushover.peakShear/results.ELF.baseShear;
    prePeakIndex = results.pushover.roofDrift < results.pushover.peakTotalDrift(nStories);
    designBaseShearDrift = interp1(results.pushover.baseShear(prePeakIndex),results.pushover.roofDrift(prePeakIndex),results.ELF.baseShear);
    effectiveYieldDrift = calcOverstr*designBaseShearDrift;
    results.pushover.periodBasedDuctility = effectiveYieldDrift/results.pushover.peak80TotalDrift(nStories);

    fprintf('Calculated overstrength = %g\n',calcOverstr)

    axis manual
    plot([0 2*max(results.pushover.roofDrift)],[results.ELF.baseShear results.ELF.baseShear],'k--')
    plot([0 results.pushover.peakTotalDrift(nStories)],[results.pushover.peakShear results.pushover.peakShear],'k--')
    plot([0 results.pushover.peak80TotalDrift(nStories)],[results.pushover.peak80Shear results.pushover.peak80Shear],'k--')

    figure
    subplot(1,2,1)
    barh(F/sum(F),0.1)
    grid on
    title('Pushover force distribution')
    ylabel('Story')
    subplot(1,2,2)
    hold on
    plot(results.pushover.peakStoryDriftRatio*100  ,1:nStories,'*-')
    plot(results.pushover.peak80StoryDriftRatio*100,1:nStories,'*-')
    grid on
    grid minor
    xl = xlim;
    xlim([0 xl(2)])
    ylim([0.5 nStories+0.5])
    ylabel('Story')
    xlabel('Story Drift Ratio (%)')
    title('Pushover story drifts')
    legend('V_{max}','V_{80}','Location','Southeast')

otherwise
    fprintf(2,'Pushover analysis failed. See results for details.\n')
    failedRatio = results.pushover.storyDrift(end,:)./bldg.storyHeight;
    figure
    subplot(1,2,1)
    barh(F/sum(F),0.1)
    grid on
    title('Pushover force distribution')
    ylabel('Story')
    subplot(1,2,2)
    hold on
    plot(failedRatio*100  ,1:nStories,'*-')
    grid on
    grid minor
    xl = xlim;
    xlim([0 xl(2)])
    ylim([0.5 nStories+0.5])
    ylabel('Story')
    xlabel('Story Drift Ratio (%)')
    title('Pushover story drifts')

end

if bldg.verbose
    pushover_time = toc(pushover_tic);
    fprintf('Pushover analysis took %g seconds.\n',pushover_time);
end

end

%##############################################################################%
%% Incremental Dynamic Analysis
if runIDA

[results.IDA,R_accepted] = incrementalDynamicAnalysis(bldg,'ground_motions.mat',results.pushover);

% Plot sample response history
randGM = randi(bldg.optionsIDA.nMotions);
randID = randi(length(bldg.optionsIDA.ST));
while isempty(results.IDA{randGM,randID})
    randGM = randi(bldg.optionsIDA.nMotions);
    randID = randi(length(bldg.optionsIDA.ST));
end
bldg.plotSampleResponse(results.IDA{randGM,randID},'story',1:nStories,'roof')

end

%##############################################################################%
%% Plot backbone curves
figure
hold on
a = zeros(nStories,1);b = zeros(nStories,1);c = zeros(nStories,1);
legendentries = cell(nStories,1);
parfor i = 1:nStories
    materialDefinition = spring(i).definition;
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';
    endpoint = spring(i).defl_y + spring(i).defl_p + spring(i).defl_pc;
    anaobj = UniaxialMaterialAnalysis(materialDefinition);
    rateType    = 'StrainRate';
    rateValue   = 0.001;
    backbone{i} = anaobj.runAnalysis([0 endpoint],rateType,rateValue,i);
    legendentries{i} = sprintf('Story %i',i);
    a(i) = spring(i).defl_y;
    b(i) = spring(i).defl_p;
    c(i) = spring(i).defl_pc;
end
results.backbone = backbone; clear backbone;
for i = 1:nStories
    plot(results.backbone{i}.disp,results.backbone{i}.force)
end
xlim([0 1.1*max(a+b+c)])
title('Backbone curves')
xlabel(sprintf('Deflection (%s)',bldg.units.length))
ylabel(sprintf('Force (%s)',bldg.units.force))
legend(legendentries)
grid on

%##############################################################################%
%% Plot hysteretic curves
if plotHysteretic
hysteretic_pos_env = cell(nStories,1);
hysteretic_neg_env = cell(nStories,1);
hysteretic = cell(nStories,1);
parfor i = 1:nStories
    materialDefinition = spring(i).definition;
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';

    anaobj = UniaxialMaterialAnalysis(materialDefinition);

    peakPoints  = [0 1 -1 2 -2 3 -3 4 -4 5 -5]*0.2*0.8*(spring(i).defl_y+spring(i).defl_p+spring(i).defl_pc);
    rateType    = 'StrainRate';
    rateValue   = peakPoints(2)/10;

    hysteretic_pos_env{i} = anaobj.runAnalysis([0 1.2*max(peakPoints)],rateType,rateValue,i);
    hysteretic_neg_env{i} = anaobj.runAnalysis([0 1.2*min(peakPoints)],rateType,rateValue,i);
    hysteretic{i}         = anaobj.runAnalysis(         peakPoints,rateType,rateValue,i);
end
results.hysteretic_pos_env = hysteretic_pos_env; clear hysteretic_pos_env;
results.hysteretic_neg_env = hysteretic_neg_env; clear hysteretic_neg_env;
results.hysteretic = hysteretic; clear hysteretic;

for i = 1:nStories
    figure
    hold on
    plot(results.hysteretic_pos_env{i}.disp,results.hysteretic_pos_env{i}.force,'k--')
    plot(results.hysteretic_neg_env{i}.disp,results.hysteretic_neg_env{i}.force,'k--')
    plot(results.hysteretic{i}.disp,results.hysteretic{i}.force,'r')
    grid on

    title(sprintf('Spring for story %i',i))
    xlabel(sprintf('Deflection (%s)',bldg.units.length))
    ylabel(sprintf('Force (%s)',bldg.units.force))
end
end
%##############################################################################%
%% Cleanup

for i = 1:length(neededPaths)
    rmpath(neededPaths{i})
end

toc
