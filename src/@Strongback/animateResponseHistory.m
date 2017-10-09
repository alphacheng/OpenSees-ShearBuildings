function animateResponseHistory(obj,results,dt)

    if nargin == 2
        dt = max(diff(results.time));
    end
    cumHeights = cumsum(obj.storyHeight);

    figure

    grid on
    grid minor

    xMax = max(max(abs(results.displacement_x)));
    yMax = sum(obj.storyHeight) + obj.storyHeight(1);

    axis([-xMax xMax 0 yMax]);
    ax = gca;
    ax.YTick = cumHeights;

    ylabel(sprintf('Height (%s)',obj.units.length))
    xlabel(sprintf('Drift (%s)',obj.units.length))

    xPos = [0 results.displacement_x(1,:)];
    yPos = [0 cumHeights] + [0 results.displacement_y(1,:)];

    h = animatedline(xPos,yPos,'Marker','*');
    t = text(2/3*xMax,1/5*yMax,'Time: 0.0s');

    for i = 1:size(results.displacement_x,1)
        displayText = sprintf('Time: %4.1fs',i*dt);
        delete(t)
        clearpoints(h)
        xPos = [0 results.displacement_x(i,:)];
        yPos = [0 cumHeights] + [0 results.displacement_y(i,:)];
        addpoints(h,xPos,yPos)
        t = text(ax,2/3*xMax,1/5*yMax,displayText);
        drawnow
    end
end
