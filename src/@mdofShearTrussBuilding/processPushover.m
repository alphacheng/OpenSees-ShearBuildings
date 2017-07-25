function results =  processPushover(obj,results,ELF)
%% PROCESSPUSHOVER Process pushover results

if ~strcmp(results.exitStatus,'Analysis Successful')
    warning('Analysis not successful; results may not be valid for post-processing')
end

results.peakShear = max(results.baseShear);
peakIndex = results.baseShear == results.peakShear;
results.peakTotalDrift = results.displacement_x(peakIndex,:);
peakStoryDrift = results.storyDrift(peakIndex,:);

postPeakIndex = results.roofDrift > results.roofDrift(peakIndex);
postPeakShear = results.baseShear(postPeakIndex);
postPeakTotalDrift = results.displacement_x(postPeakIndex,:);
postPeakStoryDrift = results.storyDrift(postPeakIndex,:);

results.peak80Shear = 0.8*results.peakShear;
results.peak80TotalDrift = interp1(postPeakShear,postPeakTotalDrift,results.peak80Shear);
peak80StoryDrift = interp1(postPeakShear,postPeakStoryDrift,results.peak80Shear);

results.peakStoryDriftRatio   = peakStoryDrift./obj.storyHeight;
results.peak80StoryDriftRatio = peak80StoryDrift./obj.storyHeight;

results.calcOverstr = results.peakShear/ELF.baseShear;
prePeakIndex = results.roofDrift < results.peakTotalDrift(obj.nStories);
designBaseShearDrift = interp1(results.baseShear(prePeakIndex),results.roofDrift(prePeakIndex),ELF.baseShear);
results.effectiveYieldDrift = results.calcOverstr*designBaseShearDrift;
results.periodBasedDuctility = results.peak80TotalDrift(obj.nStories)/results.effectiveYieldDrift;

end
