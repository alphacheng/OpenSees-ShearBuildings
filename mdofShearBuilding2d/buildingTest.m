% dynamicDrift.m

clear; close all; clc;

%% Define Building
nStories = 3;
bldg = mdofShearBuilding2d(nStories);
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;
bldg.pushover_stepSize = 0.001;
bldg.pushover_maxDrift = 6.0;

A = 90*90;          % ft^2
floorDL = 0.080;    % ksf
roofDL = 0.030;     % ksf

R = 8;
Cd = 5.5;
Omega_o = 3;

bldg.storyMass = [floorDL*A floorDL*A roofDL*A];

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
};

%% Response History Analysis
gmfile  = 'test.acc';
dt      = 0.01;
SF      = 200.0;
tend    = 50.0;
results = bldg.responseHistory(gmfile,dt,SF,tend);

% Plot sample response history
figure
subplot(211)
plot(results.time,results.groundMotion,'-')
grid on
grid minor
xlabel('Time (s)')
ylabel('Acceleration (g)')
title('Input Ground Motion')

subplot(212)
plot(results.time,results.roofDrift,'-')
grid on
grid minor
axisLimits = axis;
axis([axisLimits(1:2),-max(abs(axisLimits(3:4))),max(abs(axisLimits(3:4)))])
xlabel('Time (s)')
ylabel('Roof drift')
title('Roof Drift')

% Vary safety factor
SF = 10:10:200;
maxDrift = zeros(length(SF),1);
for i = 1:length(SF)
    results = bldg.responseHistory(gmfile,dt,SF(i),tend);
    maxDrift(i) = max(max(abs(results.storyDrift)));
end

figure
plot(maxDrift,SF,'o-')
grid on
xlabel('Maximum story drift')
ylabel('Ground motion factor')
