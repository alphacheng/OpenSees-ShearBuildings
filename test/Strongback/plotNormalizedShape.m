figure
plot([0 desired], [0 cumsum(bldg.storyHeight)], '*-')
hold on
plot([0 shape], [0 cumsum(bldg.storyHeight)], '*-')
grid on
xlim([0 1.2])
ylabel('Height (ft)')
xlabel('Normalized displacement \delta/\delta_r')
xlim([0 1])
title('Deflected shapes')
legend('Desired (no damage concentration)', 'Actual (no strongback)')
