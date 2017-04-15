% test_ELFdesign.m
% Units: kslug, kip, ft, sec

clear; close all; clc;
% load('ground_motions.mat');

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
M_y = zeros(nStories,1);
Lambda = 8;     % Cyclic deterioration parameter
Fpr = 0.2;      % Ratio of (the force at which reloading begins)/(the force corresponding to the maximum historic deformation demand)
A_Pinch = 0.2;  % Ratio of reloading stiffness
c = 1;          % rate of deterioration
theta_p = zeros(nStories,1);  % pre-capping rotation
theta_pc = 1.5; % post-capping rotation
Res = 0.2;      % residual strength ratio
theta_u = 1.9;  % ultimate rotation capacity
D = 1.0;        % rate of cyclic deterioration
nFactor = 0;    % elastic stiffness amplification factor

py_factor = 3;  % ratio of theta_p to theta_y
for i = 1:nStories
    Solution = [-py_factor 1 0 ; 0 as*K0(i) 1 ; 1 0 -1/K0(i)]^-1 * [0 ; resultsELF.designStrength(i) ; 0];
    % theta_y(i) = Solution(1); % rotation at yield
    theta_p(i) = Solution(2); % pre-capping rotation
    M_y(i)     = Solution(3); % effective yield strength
end

bldg.storySpringDefinition = cell(nStories,1);
for i = 1:nStories
    bldg.storySpringDefinition{i} = bilinearMaterialDefinition(i,K0(i),as,M_y(i),Lambda,c,theta_p(i),theta_pc,Res,theta_u,D,nFactor);
end

%% Response History Analysis
gmfile  = 'test.acc';
dt      = 0.01;
SF      = 10.0;
tend    = 33.0;
results = bldg.responseHistory(gmfile,dt,SF*bldg.g,tend,1);

% Plot sample response history
figure
subplot(211)
plot(results.time,results.groundMotion,'-')
grid on
grid minor
xlabel('Time (s)')
ylabel('Acceleration (ft/s^2)')
title('Input Ground Motion')

subplot(212)
plot(results.time,results.roofDrift,'-')
grid on
grid minor
axisLimits = axis;
axis([axisLimits(1:2),-max(abs(axisLimits(3:4))),max(abs(axisLimits(3:4)))])
xlabel('Time (s)')
ylabel('Drift (ft)')
title('Roof Drift')

% Vary safety factor
SF = 1:0.2:10;
maxDrift = zeros(length(SF),1);
parfor i = 1:length(SF)
    fprintf('Calculating IDA point for SF = %g\n',SF(i));
    results = bldg.responseHistory(gmfile,dt,SF(i)*bldg.g,tend,i);
    maxDrift(i) = max(max(abs(results.storyDrift)));
end

figure
plot(maxDrift,SF,'o-')
grid on
xlabel('Maximum story drift')
ylabel('Ground motion factor')
