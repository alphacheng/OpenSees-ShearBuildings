function fig = plotPushoverCurve(obj,results)

fig = figure;
hold on
plot(results.roofDrift,results.baseShear,'-')
grid on
grid minor
xlabel(sprintf('Roof drift (%s)',obj.units.length))
ylabel(sprintf('Base shear (%s)',obj.units.force))
title('Pushover analysis')

if nargout == 0
    clear fig
end

end
