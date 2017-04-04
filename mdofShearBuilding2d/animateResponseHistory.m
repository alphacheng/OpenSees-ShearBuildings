%% Animate response history

function animateResponseHistory( results  )

figure

grid on
grid minor

dt = max(diff(results.time));
nStories = size(results.totalDrift,2);
xMax = max(max(abs(results.totalDrift)));
yMax = nStories+1;

axis([-xMax xMax 0 yMax]);
ax = gca;
ax.YTick = 0:nStories;

ylabel('Story')
xlabel('Drift')

h = animatedline([0 results.totalDrift(1,:)],0:nStories,'Marker','*');
t = text(2/3*xMax,1/5*yMax,'Time: 0.0s');

for i = 1:size(results.totalDrift,1)
    displayText = sprintf('Time: %4.1fs',i*dt);
    delete(t)
    clearpoints(h)
    addpoints(h,[0 results.totalDrift(i,:)],0:nStories)
    t = text(ax,2/3*xMax,1/5*yMax,displayText);
    drawnow
end

end