function results =  processPushover(obj,results,ELF)
%% PROCESSPUSHOVER Process pushover results

if ~strcmp(results.exitStatus,'Analysis Successful')
    warning('Analysis not successful; results may not be valid for post-processing')
end

baseShear = results.baseShear;
roofDrift = results.roofDrift;

peakShear = max(baseShear);
peakIndex = find(baseShear == peakShear, 1);
peakTotalDrift = results.displacement_x(peakIndex,:);
peakStoryDrift = results.storyDrift(peakIndex,:);

postPeakIndex      = roofDrift > roofDrift(peakIndex);
[postPeakShear, IA] = unique(baseShear(postPeakIndex), 'stable');
postPeakTotalDrift = results.displacement_x(postPeakIndex,:);
postPeakTotalDrift = postPeakTotalDrift(IA,:);
postPeakStoryDrift = results.storyDrift(postPeakIndex,:);
postPeakStoryDrift = postPeakStoryDrift(IA,:);

peak80Shear      = 0.8*peakShear;
peak80TotalDrift = interp1(postPeakShear,postPeakTotalDrift,peak80Shear);
peak80StoryDrift = interp1(postPeakShear,postPeakStoryDrift,peak80Shear);

peakStoryDriftRatio   = peakStoryDrift./obj.storyHeight;
peak80StoryDriftRatio = peak80StoryDrift./obj.storyHeight;

prePeakIndex         = roofDrift < peakTotalDrift(obj.nStories);
[prePeakShear, IA]   = unique(baseShear(prePeakIndex), 'stable');
prePeakDrift         = roofDrift(prePeakIndex);
prePeakDrift         = prePeakDrift(IA);

calcOverstr          = peakShear/ELF.baseShear;
designBaseShearDrift = interp1(prePeakShear,prePeakDrift,ELF.baseShear);
effectiveYieldDrift  = calcOverstr*designBaseShearDrift;
periodBasedDuctility = peak80TotalDrift(obj.nStories)/effectiveYieldDrift;

% Return results
results.peakShear               = peakShear;
results.peakTotalDrift          = peakTotalDrift;
results.peakStoryDriftRatio     = peakStoryDriftRatio;
results.peak80Shear             = peak80Shear;
results.peak80TotalDrift        = peak80TotalDrift;
results.peak80StoryDriftRatio   = peak80StoryDriftRatio;
results.calcOverstr             = calcOverstr;
results.effectiveYieldDrift     = effectiveYieldDrift;
results.periodBasedDuctility    = periodBasedDuctility;

end
