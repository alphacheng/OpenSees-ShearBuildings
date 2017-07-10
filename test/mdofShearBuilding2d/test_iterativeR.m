% test_iterativeR.m

clear all; close all; clc;   %#ok<CLALL>

if isempty(gcp)
    parpool;
end

archList = [3 6 9];
tic

test_options

complete = false;
lastFalse = options.maxR;
lastTrue  = options.minR;
while complete == false

    results = struct;
    bldg = cell(length(archList),1);

    for archIndex = 1:length(archList)
        nStories = archList(archIndex);

        bldg{archIndex} = generateArchetype(nStories,options);

    %-------------------------------------------------------------------------------
    % ELF and spring definition
        results(archIndex).ELF = bldg{archIndex}.ELFanalysis();
        spring = bldg{archIndex}.springDesign(results(archIndex).ELF,springGivens);

        storySpringDefinition = cell(nStories,1);
        for i = 1:nStories
            storySpringDefinition{i} = spring(i).definition;
        end
        bldg{archIndex}.storySpringDefinition = storySpringDefinition;

    %-------------------------------------------------------------------------------
    % Pushover
        [~,eigenvecs] = bldg{archIndex}.eigenvalues;

        F = bldg{archIndex}.storyMass' .* eigenvecs;
        if min(F) < 0
            F = -F;
        end

        results(archIndex).pushover = bldg{archIndex}.pushover(F,'TargetPostPeakRatio',0.75);

        calcOverstr = results(archIndex).pushover.peakShear/results(archIndex).ELF.baseShear;
        prePeakIndex = results(archIndex).pushover.roofDrift < results(archIndex).pushover.peakTotalDrift(nStories);
        designBaseShearDrift = interp1(results(archIndex).pushover.baseShear(prePeakIndex),results(archIndex).pushover.roofDrift(prePeakIndex),results(archIndex).ELF.baseShear);
        effectiveYieldDrift = calcOverstr*designBaseShearDrift;
        results(archIndex).pushover.periodBasedDuctility = results(archIndex).pushover.peak80TotalDrift(nStories)/effectiveYieldDrift;

    %-------------------------------------------------------------------------------
    % IDA
        results(archIndex).IDA = incrementalDynamicAnalysis(bldg{archIndex},'ground_motions.mat',results(archIndex).pushover);
    end
    IDA = [results.IDA];
    ACMRavg = mean([IDA.ACMR]);
    beta_avg = mean([IDA.beta_total]);
    arch_pass = [IDA.R_accepted];
    clear IDA;

    ACMR10 = FEMAP695_ACMRxx(beta_avg,0.1);

    if ACMRavg >= ACMR10
        suite_passes = true;
    else
        suite_passes = false;
    end

% Figure out what direction to iterate
    success = ~any(~[arch_pass suite_passes]);
    if success
        lastTrue = options.respModCoeff;
        if options.respModCoeff == options.maxR
            complete = true;
        elseif (options.respModCoeff - prevR) < options.R_tol
            complete = true;
        end
    else
        lastFalse = options.respModCoeff;
        if options.respModCoeff == options.minR
            error('R has bottomed out.')
        end
    end
    if complete == false
        prevR = options.respModCoeff;
        options.respModCoeff = max(min(lastTrue + (lastFalse - lastTrue)/2, options.maxR), options.minR);
        options.deflAmplFact = options.respModCoeff;
        fprintf('\nBeginning new iteration. R = %g\n\n',options.respModCoeff)
    else
        fprintf('\nIteration complete. R = %g\n',options.respModCoeff)
    end
    % check = [arch_pass suite_passes];
    % if any(~check)
    %     options.respModCoeff = options.respModCoeff - 0.5;
    %     options.deflAmplFact = options.deflAmplFact - 0.5;
    % else
    %     if R == 8
    %         complete = true;
    %     else
    %         alreadyBeenHere = any()
    %         if alreadyBeenHere
    %             complete = true;
    %         else
    %             options.respModCoeff = options.respModCoeff + 0.5;
    %             options.deflAmplFact = options.deflAmplFact + 0.5;
    %         end
    %     end
    % end
end

toc
