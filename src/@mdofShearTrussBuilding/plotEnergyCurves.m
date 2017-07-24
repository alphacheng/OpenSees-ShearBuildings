function fig = plotEnergyCurves(obj,rh)
fig = figure;
plot(rh.time,[rh.energy.earthquake rh.energy.norm_gravity])
yl = ylim;
ylim([0 yl(2)])
grid on
xlabel(sprintf('Time (%s)',obj.units.time))
ylabel(sprintf('Energy (%s*%s)',obj.units.length,obj.units.force))
if nargout == 0
    clear fig
end
end
