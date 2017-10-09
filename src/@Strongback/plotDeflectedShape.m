function plotDeflectedShape(obj, results, index)

cumHeights = cumsum(obj.storyHeight);

figure

xPos = [0 results.displacement_x(index,:)];
yPos = [0 cumHeights] + [0 results.displacement_y(index,:)];

plot(xPos, yPos, 'k*-')

grid on
grid minor

xMax = max(max(abs(results.displacement_x)));
yMax = sum(obj.storyHeight) + obj.storyHeight(1);

axis([-xMax xMax 0 yMax]);
ax = gca;
ax.YTick = cumHeights;

ylabel(sprintf('Height (%s)',obj.units.length))
xlabel(sprintf('Drift (%s)',obj.units.length))

end
