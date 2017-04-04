% pushoverDrift.m

clear; close all; clc;

%% Define Building
nStories = 4;
bldg = mdofShearBuilding2d(nStories);
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;
bldg.pushover_stepSize = 0.001;
bldg.pushover_maxDrift = 6.0;

bldg.storyMass = [2 1 1 1];

% Define Material
K0 = 1000;      % elastic stiffness
as = 0.04;      % strain hardening ratio
My = 100;       % effective yield strength
Lambda = 8;     % Cyclic deterioration parameter
Fpr = 0.2;      % Ratio of (the force at which reloading begins)/(the force corresponding to the maximum historic deformation demand)
A_Pinch = 0.2;  % Ratio of reloading stiffness
c = 1;          % rate of deterioration
theta_p = 0.3;  % pre-capping rotation
theta_pc = 1.5; % post-capping rotation
Res = 0.2;      % residual strength ratio
theta_u = 1.9;  % ultimate rotation capacity
D = 1.0;        % rate of cyclic deterioration
nFactor = 0;    % elastic stiffness amplification factor

materialDefinition = bilinearMaterialDefinition(K0,as,My,Lambda,c,theta_p,theta_pc,Res,theta_u,D,nFactor);
% materialDefinition = peakOrientedMaterialDefinition(K0,as,My,Lambda,c,theta_p,theta_pc,Res,theta_u,D);
% materialDefinition = pinchingMaterialDefinition(K0,as,My,Fpr,A_Pinch,Lambda,c,theta_p,theta_pc,Res,theta_u,D);

bldg.storySpringDefinition = {
    materialDefinition
    sprintf('uniaxialMaterial Elastic 2 %g',K0)
    sprintf('uniaxialMaterial Elastic 3 %g',K0)
    sprintf('uniaxialMaterial Elastic 4 %g',K0)
};

%% Pushover Analysis
results = bldg.pushover(10*[0.1 0.3 0.3 0.3]','TargetPostPeakRatio',0.75);

peakShear = max(results.baseShear);
peakShearIndex = results.baseShear == peakShear;

postPeakIndex = results.roofDrift > results.roofDrift(peakShearIndex);

peakShear80 = 0.8*peakShear;
peakShear80Index = find(postPeakIndex & (results.baseShear < peakShear80),1);

figure
hold on
plot(results.roofDrift,results.baseShear,'-')
grid on
grid minor
axis manual
ax = gca;
noteText = sprintf('80%% Peak Shear = %g',peakShear80);
text(results.roofDrift(peakShear80Index),1.1*peakShear80,noteText)
plot([0 ax.XLim(2)],[peakShear80 peakShear80],'k-')
xlabel('Roof drift')
ylabel('Base shear')
title('P-Delta for Roof')

%% Deflected shape - Pushover
figure
% Peak Shear
subplot(1,2,1)
plot([0 results.totalDrift(peakShearIndex,:)],0:nStories,'*-')
grid on
grid minor
ylabel('Story')
xlabel('Total drift')
title('Peak Shear')

% 80% of Peak Shear
subplot(1,2,2)
plot([0 results.totalDrift(peakShear80Index,:)],0:nStories,'*-')
grid on
grid minor
ylabel('Story')
xlabel('Total drift')
title('80% of Peak Shear')
