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

results = struct;

%##############################################################################%
%% Define Building
nStories = 3;
bldg = mdofShearBuilding2d(nStories);

% Units
bldg.units.force = 'kip';
bldg.units.mass  = 'kslug';
bldg.units.length= 'ft';
bldg.units.time  = 'sec';

bldg.g = 32.2;                                  % Acceleration due to gravity

bldg.storyHeight = [20 15 15];                  % Story heights (ft)
storyDL = [80 80 30]/1000;                  % Story dead loads (ksf)
storyArea = 90*90*ones(1,nStories);             % Story areas (ft^2)
bldg.storyMass = (storyDL .* storyArea)/bldg.g; % Story masses (kslug)

bldg.seismicDesignCategory = 'Dmax';
bldg.respModCoeff = 8;
bldg.deflAmplFact = 8;
bldg.overstrengthFactor = 3;
bldg.impFactor = 1;

springGivens.as       =  0.05;  % strain hardening ratio
springGivens.Lambda_S =  8.00;  % Cyclic deterioration parameter - strength
springGivens.Lambda_K = 10.00;  % Cyclic deterioration parameter - stiffness
springGivens.c_S      =  1.00;  % rate of deterioration - strength
springGivens.c_K      =  1.00;  % rate of deterioration - stiffness
springGivens.Res      =  0.30;  % residual strength ratio
springGivens.D        =  1.00;  % rate of cyclic deterioration
springGivens.nFactor  =  0.00;  % elastic stiffness amplification factor
springGivens.C_yc     =  0.80;  % ratio of yield strength to capping strength
springGivens.C_pcp    =  1.00;  % ratio of post-capping deflection to pre-capping deflection
springGivens.C_upc    = 20.00;  % ratio of ultimate deflection to u_y + u_p + u_pc

springGivens.stiffnessSafety = 1.0;

springGivens.enforceMinimum = true;
springGivens.minimumRatio = 0.6;


%##############################################################################%
%% Analysis Options
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = true;

verbose     = true ;    % Toggle verbose output
runPushover = true ;    % Toggle pushover analysis
runIDA      = false;    % Toggle IDA

% Equivalent lateral force options
iterate = false;             % Select whether to do iteration
iterOption = 'overstrength';      % Variable to use for convergence: 'period' or 'overstrength'

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
nIter = 0;
while iterating == true
    nIter = nIter + 1;
    %% Define springs
    results.ELF = bldg.ELFanalysis();
    spring = bldg.springDesign(results.ELF,springGivens);

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
                    fprintf('Did %i iterations.\n',nIter)
                end
            case 'overstrength'
                % [~,eigenvecs] = bldg.eigenvalues;
                [~,modeShapes] = modalAnalysis(bldg,spring);
                eigenvecs = modeShapes(1,:)';
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

% [~,eigenvecs] = bldg.eigenvalues;
[~,modeShapes] = modalAnalysis(bldg,spring);
eigenvecs = modeShapes(1,:)';

F = bldg.storyMass' .* eigenvecs;
if min(F) < 0
    error('Mode shapes invalid.')
end

results.pushover = bldg.pushover(F,'TargetPostPeakRatio',0.75);

switch results.pushover.exitStatus
case 'Analysis Successful'
    figure
    hold on
    plot(results.pushover.roofDrift,results.pushover.baseShear,'-')
    axis manual
    plot([0 2*max(results.pushover.roofDrift)],[results.ELF.baseShear results.ELF.baseShear],'k--')
    plot([0 results.pushover.peakTotalDrift(nStories)],[results.pushover.peakShear results.pushover.peakShear],'k--')
    plot([0 results.pushover.peak80TotalDrift(nStories)],[results.pushover.peak80Shear results.pushover.peak80Shear],'k--')
    grid on
    grid minor
    xlabel(sprintf('Roof drift (%s)',bldg.units.length))
    ylabel(sprintf('Base shear (%s)',bldg.units.force))
    title('Pushover analysis')

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
    ylim([0.5 nStories+0.5])
    ylabel('Story')
    xlabel('Story Drift Ratio (%)')
    title('Pushover story drifts')
    legend('V_{max}','V_{80}','Location','Southeast')

otherwise
    fprintf(2,'Pushover analysis failed. See results for details.\n')

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
        maxDriftRatio(j) = max(max(abs(IDA{i,j}.storyDrift))./bldg.storyHeight);
        if verbose
            fprintf('Maximum story drift ratio = %5.2f%%\n',maxDriftRatio(j)*100);
        end
    end

    plot(maxDriftRatio*100,ST,'o-')
    legendentries{i} = ground_motions(i).ID;

    if bldg.deleteFilesAfterAnalysis
        delete(gmfile)
    end
end
results.IDA = IDA;

grid on
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
% backbone_ax = cell(nStories,1);
legendentries = cell(nStories,1);
for i = 1:nStories
    % backbone_ax{i} = subplot(nStories,1,i);
    materialDefinition = bldg.storySpringDefinition{i};
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';
    plotBackboneCurve(materialDefinition,spring.defl_u(i),false)
    legendentries{i} = sprintf('Story %i',i);
end
xlim([0 1.1*max(spring.defl_y+spring.defl_p+spring.defl_pc)])
title('Backbone curves')
xlabel(sprintf('Deflection (%s)',bldg.units.length))
ylabel(sprintf('Force (%s)',bldg.units.force))
legend(legendentries)
grid on

%##############################################################################%
%% Plot hysteretic curves
for i = 1:nStories
    materialDefinition = bldg.storySpringDefinition{i};
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';

    anaobj = UniaxialMaterialAnalysis(materialDefinition);

    peakPoints  = [0 1 -1 2 -2 3 -3 4 -4 5 -5]*0.2;
    rateType    = 'StrainRate';
    rateValue   = peakPoints(2)/10;

    results.hysteretic_pos_env = anaobj.runAnalysis([0 max(peakPoints)],rateType,rateValue);
    results.hysteretic_neg_env = anaobj.runAnalysis([0 min(peakPoints)],rateType,rateValue);
    results.hysteretic         = anaobj.runAnalysis(         peakPoints,rateType,rateValue);

    figure
    hold on
    plot(results.hysteretic_pos_env.disp,results.hysteretic_pos_env.force,'k--')
    plot(results.hysteretic_neg_env.disp,results.hysteretic_neg_env.force,'k--')
    plot(results.hysteretic.disp,results.hysteretic.force,'r')
    grid on
end

%##############################################################################%
%% Cleanup

for i = 1:length(neededPaths)
    rmpath(neededPaths{i})
end

toc
