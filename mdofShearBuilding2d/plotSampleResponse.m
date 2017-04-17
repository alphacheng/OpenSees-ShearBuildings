function plotSampleResponse(results)
% plotSampleResponse Description
%    plotSampleResponse(results)
%
% Long description
%
%

figure
subplot(211)
plot(results.time,results.groundMotion,'-')
grid on
grid minor
xlabel('Time (s)')
ylabel('Acceleration (ft/s^2)')
titleText = sprintf('Input Ground Motion (GM: %s, Index: %i)',results.gmID,results.indexNum);
title(titleText)

subplot(212)
plot(results.time,results.roofDrift,'-')
grid on
grid minor
axisLimits = axis;
axis([axisLimits(1:2),-max(abs(axisLimits(3:4))),max(abs(axisLimits(3:4)))])
xlabel('Time (s)')
ylabel('Drift (ft)')
title('Roof Drift')


% function end: 'plotSampleResponse'
end
