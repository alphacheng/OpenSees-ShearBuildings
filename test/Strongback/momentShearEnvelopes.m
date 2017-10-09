function h = momentShearEnvelopes(bldg, analysisType, ground_motions, gmIndex)
%MOMENTSHEARENVELOPES  Calculate and plot moment and shear envelopes
%
%   h = momentShearEnvelopes(obj, analysisType, ground_motions, gmIndex)
%
%   ground_motions and gmIndex are only required if 'responseHistory' is
%   selected as the analysisType.
%

% Create moment plot figure
h(1) = figure;
ax_m = gca;
hold on
grid on
% Create shear plot figure
h(2) = figure;
ax_s = gca;
hold on
grid on
switch analysisType
case 'responseHistory'
    for gmi = gmIndex(1):gmIndex(2)
        gmID   = ground_motions(gmi).ID;
        gmFile = scratchFile(bldg,sprintf('acc%s.acc',gmID));
        dt     = ground_motions(gmi).dt;
        tEnd   = ground_motions(gmi).time(end) + bldg.optionsIDA.tExtra;
        SF = 1;
        accel = ground_motions(gmi).normalized_acceleration;
        dlmwrite(gmFile,accel*bldg.g);

        results = bldg.responseHistory(gmFile, dt, SF, tEnd, gmID, gmi);

        % Plot moment envelope
        maxMoment = max(abs(results.strongbackMoment));
        plot(ax_m, maxMoment, [0 cumsum(bldg.storyHeight)], '*-')

        % Plot shear envelope
        maxShear = max(abs(results.strongbackShear));
        plot(ax_s, maxShear, [0 cumsum(bldg.storyHeight)], '*-')
    end
    title(ax_m, sprintf('Moment envelope, %i-story building\n EI=%g %s-%s^2, gm %s-%s',...
        bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia, bldg.units.force,...
        bldg.units.length, ground_motions(gmIndex(1)).ID, ground_motions(gmIndex(2)).ID))
    title(ax_s, sprintf('Shear envelope, %i-story building\n EI=%g %s-%s^2, gm %s-%s',...
        bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia, bldg.units.force,...
        bldg.units.length, ground_motions(gmIndex(1)).ID, ground_motions(gmIndex(2)).ID))
    legend(ax_m, {ground_motions(gmIndex(1):gmIndex(2)).ID})
    legend(ax_s, {ground_motions(gmIndex(1):gmIndex(2)).ID})
case 'pushover'
    F = bldg.pushoverForceDistribution();
    results = bldg.pushover(F, 'TargetPostPeakRatio', 0.75);

    % Plot moment envelope
    maxMoment = max(abs(results.strongbackMoment));
    plot(ax_m, maxMoment, [0 cumsum(bldg.storyHeight)], '*-')

    % Plot shear envelope
    maxShear = max(abs(results.strongbackShear));
    plot(ax_s, maxShear, [0 cumsum(bldg.storyHeight)], '*-')

    title(ax_m, sprintf('Moment envelope, %i-story building\n EI=%g %s-%s^2, pushover',...
        bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia,...
        bldg.units.force, bldg.units.length))
    title(ax_s, sprintf('Shear envelope, %i-story building\n EI=%g %s-%s^2, pushover',...
        bldg.nStories, bldg.strongbackDefinition.Modulus*bldg.strongbackDefinition.Inertia,...
        bldg.units.force, bldg.units.length))
end
xlabel(ax_m, sprintf('Moment (%s-%s)', bldg.units.force, bldg.units.length))
xlabel(ax_s, sprintf('Shear (%s)', bldg.units.force))
ylabel(ax_m, sprintf('Height (%s)', bldg.units.length))
ylabel(ax_s, sprintf('Height (%s)', bldg.units.length))

end
