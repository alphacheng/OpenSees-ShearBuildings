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

%##############################################################################%
%% Equivalent lateral force procedure

if verbose; fprintf('Running equivalent lateral force procedure...\n'); end

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
if verbose
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
IDA = cell(nMotions,length(SF2));
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
            fprintf('Calculating IDA point for {%s, %2i}, SF2 = %5.2f ... ',ground_motions(i).ID,j,SF2(j));
        end
        SF = SF1*SF2(j);
        IDA{i,j} = bldg.responseHistory(gmfile,dt,SF,tend,ground_motions(i).ID,j);
        switch IDA{i,j}.exitStatus
        case 'Analysis Failed'
            maxDriftRatio(j) = NaN;
            if verbose
                fprintf('Analysis failed\n');
            end
        case 'Analysis Successful'
            maxDriftRatio(j) = max(max(abs(IDA{i,j}.storyDrift))./bldg.storyHeight);
            if verbose
                fprintf('Maximum story drift ratio = %5.2f%%\n',maxDriftRatio(j)*100);
            end
        end

    end
    goodDrifts = ~isnan(maxDriftRatio);

    plot(maxDriftRatio(goodDrifts)*100,ST(goodDrifts),'o-')
    legendentries{i} = ground_motions(i).ID;

    if bldg.deleteFilesAfterAnalysis
        delete(gmfile)
    end
end
results.IDA = IDA;

grid on
xlim([0 15])
xlabel('Maximum story drift ratio (%)')
ylabel('Ground motion intensity, S_T (g)')
legend(legendentries)

% Plot sample response history
plotSampleResponse(results.IDA{1,5})

if verbose
    ida_time = toc(ida_tic);
    fprintf('Incremental dynamic analysis took %g seconds.\n',ida_time);
end

end

%##############################################################################%
%% Plot backbone curves
figure
hold on
a = zeros(nStories,1);b = zeros(nStories,1);c = zeros(nStories,1);
legendentries = cell(nStories,1);
for i = 1:nStories
    materialDefinition = bldg.storySpringDefinition{i};
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';
    plotBackboneCurve(materialDefinition,spring(i).defl_u,false)
    legendentries{i} = sprintf('Story %i',i);
    a(i) = spring(i).defl_y;
    b(i) = spring(i).defl_p;
    c(i) = spring(i).defl_pc;
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
for i = 1:nStories
    materialDefinition = bldg.storySpringDefinition{i};
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';

    anaobj = UniaxialMaterialAnalysis(materialDefinition);

    peakPoints  = [0 1 -1 2 -2 3 -3 4 -4 5 -5]*0.2*0.8*(spring(i).defl_y+spring(i).defl_p+spring(i).defl_pc);
    rateType    = 'StrainRate';
    rateValue   = peakPoints(2)/10;

    results.hysteretic_pos_env = anaobj.runAnalysis([0 1.2*max(peakPoints)],rateType,rateValue);
    results.hysteretic_neg_env = anaobj.runAnalysis([0 1.2*min(peakPoints)],rateType,rateValue);
    results.hysteretic         = anaobj.runAnalysis(         peakPoints,rateType,rateValue);

    figure
    hold on
    plot(results.hysteretic_pos_env.disp,results.hysteretic_pos_env.force,'k--')
    plot(results.hysteretic_neg_env.disp,results.hysteretic_neg_env.force,'k--')
    plot(results.hysteretic.disp,results.hysteretic.force,'r')
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
