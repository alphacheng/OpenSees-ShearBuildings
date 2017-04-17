%% Animate response history

function animateResponseHistory( results, heights  )

cumHeights = zeros(1,length(heights));
for i = 1:length(heights)
    cumHeights(i) = sum(heights(1:i));
end

figure

grid on
grid minor

dt = max(diff(results.time));
xMax = max(max(abs(results.totalDrift)));
yMax = sum(heights) + heights(1);

axis([-xMax xMax 0 yMax]);
ax = gca;
ax.YTick = cumHeights;

ylabel('Height')
xlabel('Drift')

h = animatedline([0 results.totalDrift(1,:)],[0 cumHeights],'Marker','*');
t = text(2/3*xMax,1/5*yMax,'Time: 0.0s');

for i = 1:size(results.totalDrift,1)
    displayText = sprintf('Time: %4.1fs',i*dt);
    delete(t)
    clearpoints(h)
    addpoints(h,[0 results.totalDrift(i,:)],[0 cumHeights])
    t = text(ax,2/3*xMax,1/5*yMax,displayText);
    drawnow
end

end
