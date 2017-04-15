% test_ELFdesign.m
% Units: kslug, kip, ft, sec

clear; close all; clc;
load('ground_motions.mat');

%% Define Building
nStories = 3;
bldg = mdofShearBuilding2d(nStories);
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = false;
bldg.g = 32.2;

bldg.storyHeight = [20 15 15];
storyDL = [0.080 0.080 0.030];
storyArea = 90*90*ones(1,nStories);
bldg.storyMass = (storyDL .* storyArea)/bldg.g;

%% Design building
bldg.seismicDesignCategory = 'Dmax';
bldg.responseModificationCoefficient = 8;
bldg.deflectionAmplificationFactor = 5.5;
bldg.overstrengthFactor = 3;
bldg.importanceFactor = 1;

resultsELF = bldg.ELFdesign();

%% Define springs
K0 = resultsELF.designStiffness;      % elastic stiffness
as = 0.04;      % strain hardening ratio
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

[theta_y theta_p M_y] = [-py_factor 1 0 ; 0 as*K0 1 ; 1 0 -1/K0]^-1 * [0 ; resultsELF.designStrength ; 0]

My = resultsELF.designStrength - as*K0*theta_p;       % effective yield strength

bldg.storySpringDefinition = cell(nStories,1);
for i = 1:nStories
    bldg.storySpringDefinition{i} = bilinearMaterialDefinition(K0(i),as,My(i),Lambda,c,theta_p,theta_pc,Res,theta_u,D,nFactor);
end

%% Response History Analysis
gmfile  = 'test.acc';
dt      = 0.01;
SF      = 200.0;
tend    = 33.0;
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
