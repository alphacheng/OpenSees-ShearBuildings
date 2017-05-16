%##############################################################################%
%% Define Building
nStories = 3;
bldg = mdofShearBuilding2d(nStories);

% Units
bldg.units.force = 'kip';
bldg.units.mass  = 'kslug';
bldg.units.length= 'ft';
bldg.units.time  = 'sec';

bldg.g = 32.2;                                  % Acceleration due to gravity

bldg.storyHeight = [20 15 15];                  % Story heights (ft)
storyDL = [80 80 30]/1000;                  % Story dead loads (ksf)
storyArea = 90*90*ones(1,nStories);             % Story areas (ft^2)
bldg.storyMass = (storyDL .* storyArea)/bldg.g; % Story masses (kslug)

bldg.seismicDesignCategory = 'Dmax';
bldg.respModCoeff = 8;
bldg.deflAmplFact = 5;
bldg.overstrengthFactor = 3;
bldg.impFactor = 1;

springGivens.as       =  0.05;  % strain hardening ratio
springGivens.Lambda_S =  8.00;  % Cyclic deterioration parameter - strength
springGivens.Lambda_K = 10.00;  % Cyclic deterioration parameter - stiffness
springGivens.c_S      =  1.00;  % rate of deterioration - strength
springGivens.c_K      =  1.00;  % rate of deterioration - stiffness
springGivens.Res      =  0.30;  % residual strength ratio
springGivens.D        =  1.00;  % rate of cyclic deterioration
springGivens.nFactor  =  0.00;  % elastic stiffness amplification factor
springGivens.C_yc     =  0.80;  % ratio of yield strength to capping strength
springGivens.C_pcp    =  1.00;  % ratio of post-capping deflection to pre-capping deflection
springGivens.C_upc    = 20.00;  % ratio of ultimate deflection to u_y + u_p + u_pc

springGivens.stiffnessSafety = 1.0;
springGivens.strengthSafety = 1.0;

springGivens.enforceMinimum = true;
springGivens.minimumRatio = 0.5;


%##############################################################################%
%% Analysis Options
bldg.echoOpenSeesOutput = false;
bldg.deleteFilesAfterAnalysis = true;

verbose     = true ;    % Toggle verbose output
runPushover = true ;    % Toggle pushover analysis
runIDA      = true ;    % Toggle IDA

% Equivalent lateral force options
iterate = false;             % Select whether to do iteration
iterOption = 'overstrength';      % Variable to use for convergence: 'period' or 'overstrength'

diffTol = 1e-3;             % Tolerance for iteration

% Pushover analysis options
bldg.pushover_stepSize   = 0.001;
bldg.pushover_maxDrift   = 100;

% Incremental dynamic analysis options
nMotions = 12;                              % Number of ground motions to analyze
SF2 = [0:0.25:1.5 , 2:0.5:5 , 5.75:0.75:8]; % Scale factors to use for each IDA curve
