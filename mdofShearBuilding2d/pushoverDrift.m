% pushoverDrift.m

clear; close all; clc;

%% Define Building
nStories = 4;
bldg = mdofShearBuilding2d(nStories);
bldg.echoOpenSeesOutput = true;
bldg.deleteFilesAfterAnalysis = false;
bldg.pushover_stepSize = 0.001;
bldg.pushover_maxDrift = 6.0;

bldg.storyMass = [2 1 1 1];

% Define Bilinear Material
K0 = 1000;      % elastic stiffness
as = 0.04;      % strain hardening ratio
My = 100;       % effective yield strength
Lambda = 8;     % Cyclic deterioration parameter
c = 1;          % rate of deterioration
theta_p = 0.3;  % pre-capping rotation
theta_pc = 1.5; % post-capping rotation
Res = 0.2;      % residual strength ratio
theta_u = 1.9;  % ultimate rotation capacity
D = 1.0;        % rate of cyclic deterioration
nFactor = 0;    % elastic stiffness amplification factor

materialDefinition = bilinearMaterialDefinition(K0,as,My,Lambda,c,theta_p,theta_pc,Res,theta_u,D,nFactor);

bldg.storySpringDefinition = {
    materialDefinition
    'uniaxialMaterial Elastic 2 1000'
    'uniaxialMaterial Elastic 3 1000'
    'uniaxialMaterial Elastic 4 1000'
};

%% Pushover Analysis
results = bldg.pushover(10*[0.1 0.3 0.3 0.3]','TargetPostPeakRatio',0.75);

figure
plot(results1.roofDrift,results1.baseShear,'-')
grid on
grid minor
xlabel('Roof drift')
ylabel('Base shear')
title('P-Delta for Roof')

%% Deflected shape - Pushover
figure
plot([0 results1.totalDrift(end,:)],0:nStories,'*-')
grid on
grid minor
ylabel('Story')
xlabel('Total drift')
title('Deflected Shape - Pushover Analysis')
