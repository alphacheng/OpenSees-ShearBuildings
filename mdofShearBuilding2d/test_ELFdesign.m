% test_ELFdesign.m
% Units: kslug, kip, ft, sec

tic
clear; close all; clc;
addpath('../UniaxialMaterialAnalysis');

%% Define Building
nStories = 3;
bldg = mdofShearBuilding2d(nStories);
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = true;
verbose = true;                                 % Verbose output (bool)
bldg.g = 32.2;                                  % Acceleration due to gravity (ft/s^2)

bldg.storyHeight = [20 15 15];                  % Story heights (ft)
storyDL = [0.080 0.080 0.030];                  % Story dead loads (ksf)
storyArea = 90*90*ones(1,nStories);             % Story areas (ft^2)
bldg.storyMass = (storyDL .* storyArea)/bldg.g; % Story masses (kslug)

%% Design building
bldg.seismicDesignCategory = 'Dmax';
bldg.responseModificationCoefficient = 8;
bldg.deflectionAmplificationFactor = 5.5;
bldg.overstrengthFactor = 3;
bldg.importanceFactor = 1;

resultsELF = bldg.ELFdesign();

periodDiff = 1;
minDiff = 1e-3;
while abs(periodDiff) > minDiff
    %% Define springs
    K0 = resultsELF.designStiffness;    % elastic stiffness
    as = 0.04;                          % strain hardening ratio
    M_y = zeros(nStories,1);            % effective yield strength
    Lambda = 8;                         % Cyclic deterioration parameter
    c = 1;                              % rate of deterioration
    theta_p = zeros(nStories,1);        % pre-capping rotation
    Res = 0.2;                          % residual strength ratio
    D = 1.0;                            % rate of cyclic deterioration
    nFactor = 0;                        % elastic stiffness amplification factor

    py_factor = 1;  % ratio of theta_p to theta_y
    theta_y = zeros(nStories,1);
    for i = 1:nStories
        Solution = [-py_factor 1 0 ; 0 as*K0(i) 1 ; 1 0 -1/K0(i)]^-1 * [0 ; resultsELF.designStrength(i) ; 0];
        theta_y(i) = Solution(1); % rotation at yield
        theta_p(i) = Solution(2); % pre-capping rotation
        M_y(i)     = Solution(3); % effective yield strength
    end

    theta_pc = 2*theta_p;                       % post-capping rotation
    theta_u = theta_y + theta_p + 2/3*theta_pc; % ultimate rotation capacity

    bldg.storySpringDefinition = cell(nStories,1);
    for i = 1:nStories
        bldg.storySpringDefinition{i} = bilinearMaterialDefinition(i,K0(i),as,M_y(i),Lambda,c,theta_p(i),theta_pc(i),Res,theta_u(i),D,nFactor);
    end

    eigenvals = bldg.eigenvalues();
    calculatedPeriod = (eigenvals(1)/(2*pi))^-1;
    prevDiff = periodDiff;
    periodDiff = bldg.fundamentalPeriod - calculatedPeriod;
    bldg.fundamentalPeriod = calculatedPeriod;
    resultsELF = bldg.ELFdesign;
    if prevDiff == periodDiff
        break
    end
end

%% Response History Analysis
load('ground_motions.mat');
SF2 = [0:.25:1.5, 2:0.5:6, 6.75:0.75:9, 10:12];
SMT = FEMAP695_SMT(bldg.fundamentalPeriod,bldg.seismicDesignCategory);
ST  = SMT*SF2;
figure
hold on
nMotions = 6;
legendentries = cell(nMotions,1);
results = cell(nMotions,length(SF2));
for i = 1:nMotions %length(ground_motions)
    gmfile = bldg.scratchFile(sprintf('acc%s.acc',ground_motions(i).ID));
    gmfid = fopen(gmfile,'w+');
    for k = 1:ground_motions(i).numPoints
        fprintf(gmfid,'%f\n',ground_motions(i).normalized_acceleration(k)*bldg.g);
    end
    fclose(gmfid);
    dt      = ground_motions(i).dt;
    SF1     = FEMAP695_SF1(bldg.fundamentalPeriod,bldg.seismicDesignCategory);
    tend    = max(ground_motions(i).time) + 5;

    % Vary scale factor

    maxDrift = zeros(length(SF2),1);
    parfor j = 1:length(SF2)
        if verbose
            fprintf('Calculating IDA point for %s, SF2 = %g\n',ground_motions(i).ID,SF2(j));
        end
        SF = SF1*SF2(j);
        results{i,j} = bldg.responseHistory(gmfile,dt,SF,tend,ground_motions(i).ID,j);
        maxDrift(j) = max(max(abs(results{i,j}.storyDrift))./bldg.storyHeight);
    end

    plot(maxDrift,ST,'o-')
    legendentries{i} = ground_motions(i).ID;

    if bldg.deleteFilesAfterAnalysis
        delete(gmfile)
    end
end


grid on
xlabel('Maximum story drift ratio')
ylabel('Ground motion intensity, S_T (g)')
legend(legendentries)

% Plot sample response history
plotSampleResponse(results{1,5})

% Plot backbone curves
figure
hold on
for i = 1:nStories
    subplot(nStories,1,i)
    materialDefinition = bldg.storySpringDefinition{i};
    matTagLoc = strfind(materialDefinition,num2str(i));
    materialDefinition(matTagLoc(1)) = '1';
    plotBackboneCurve(materialDefinition,theta_u(i),false)
end

rmpath('../UniaxialMaterialAnalysis');

toc
