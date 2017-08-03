%##############################################################################%
%% Define Building
% nStories = 6;
% bldg = mdofShearBuilding2d(nStories);

% Units
options.units.force = 'kip';
options.units.mass  = 'kslug';
options.units.length= 'ft';
options.units.time  = 'sec';

options.g = 32.2;                                  % Acceleration due to gravity

options.storyHeight = 15;         % Story heights (ft)
options.firstHeight = 20;
options.storyDL     = 0.080;               % Story dead loads (ksf)
options.roofDL      = 0.080;
options.storyArea   = 90*90;           % Story areas (ft^2)

options.seismicDesignCategory = 'Dmax';
options.respModCoeff = 8;
options.deflAmplFact = 8;
options.overstrengthFactor = 3;
options.impFactor = 1;

springGivens.as       =  0.03;  % strain hardening ratio
springGivens.Lambda_S = 10.00;  % Cyclic deterioration parameter - strength
springGivens.Lambda_K = 10.00;  % Cyclic deterioration parameter - stiffness
springGivens.c_S      =  1.00;  % rate of deterioration - strength
springGivens.c_K      =  1.00;  % rate of deterioration - stiffness
springGivens.Res      =  0.30;  % residual strength ratio (relative to yield)
springGivens.D        =  1.00;  % rate of cyclic deterioration
springGivens.nFactor  =  0.00;  % elastic stiffness amplification factor
springGivens.C_yc     =  0.80;  % ratio of yield strength to capping strength
springGivens.C_pcp    =  2.00;  % ratio of post-capping deflection to pre-capping deflection
springGivens.C_upc    = 20.00;  % ratio of ultimate deflection to u_y + u_p + u_pc
springGivens.theta_pc = 0.3;     % angle of post-capping stiffness (degrees)
springGivens.ad       = 0.1;    % deterioration stiffness ratio
springGivens.includePDelta = false;

springGivens.stiffnessSafety = 1.0;
springGivens.strengthSafety  = 1.0;

springGivens.enforceMinimumStiffness = false;
springGivens.enforceMinimumStrength = false;
springGivens.minimumRatio = 0.7;

%##############################################################################%
%% Analysis Options
options.echoOpenSeesOutput = false;
options.deleteFilesAfterAnalysis = true;

options.verbose   = true ; % Toggle verbose output
runPushover    = true ; % Toggle pushover analysis
runIDA         = true ; % Toggle IDA
plotHysteretic = false; % Toggle plotting hysteretic curves

options.minR = 1;
options.maxR = 8;
options.Rtol = 0.5;

%-------------------------------------------------------------------------------
% Equivalent lateral force options
iterate    = false;             % Select whether to do iteration
iterOption = 'overstrength';    % Variable to use for convergence: 'period' or 'overstrength'

diffTol = 1e-3;             % Tolerance for iteration

%-------------------------------------------------------------------------------
% Pushover analysis options
options.pushover_stepSize = 0.001;
options.pushover_maxDrift = min(options.storyHeight);

options.pushover.constraints.type = 'Plain';
options.pushover.constraints.penalty.alphaS = 1.0e12;
options.pushover.constraints.penalty.alphaM = 1.0e12;

options.pushover.test.type       = 'NormDispIncr';
options.pushover.test.tolerance  = [1e-5,1e-4,1e-3];
options.pushover.test.iterations = 10;
options.pushover.test.print      = 0;
options.pushover.test.normType   = 2;

options.pushover.algorithm.type = { 'KrylovNewton','Newton','ModifiedNewton' };

%-------------------------------------------------------------------------------
% Response history options
options.damping_ModeA  = 1;             % Mode A for rayleigh damping
options.damping_ModeB  = 3;             % Mode B for rayleigh damping
options.damping_RatioA = 0.05;          % Damping ratio for mode A
options.damping_RatioB = 0.05;          % Damping ratio for mode B

options.responseHistory.constraints.type = 'Transformation';
options.responseHistory.constraints.penalty.alphaS = 1.0e12;
options.responseHistory.constraints.penalty.alphaM = 1.0e12;

options.responseHistory.test.type       = 'NormDispIncr';
options.responseHistory.test.tolerance  = [1e-5,1e-4,1e-3];
options.responseHistory.test.iterations = 10;
options.responseHistory.test.print      = 0;
options.responseHistory.test.normType   = 2;

options.responseHistory.algorithm.type = { 'KrylovNewton','Newton','ModifiedNewton' };

%-------------------------------------------------------------------------------
% Incremental dynamic analysis options
useCustomMotions = false;
RSNS = [953 960 1602 1787 169 174 1111 1116 1158 1148 900 848 752 767 1633 721 725 829 1244 1485 68 125];
spectra = '/home/petertalley/Downloads/Updated_NGA_West2_flatfiles_part1/NGA_West2_5%_Damped_Spectral_Intensity.csv';
options.IDA.nMotions = 7;                              % Number of ground motions to analyze
options.IDA.tExtra = 5;                                 % Extra analysis time after end of ground motion
options.IDA.collapseDriftRatio = 0.05;                  % Story drift ratio that defines collapse
options.IDA.ST = [0.25:0.125:4 , 4.25:0.25:8];
options.IDA.rating_DR = 'C';
options.IDA.rating_TD = 'C';
options.IDA.rating_MDL = 'C';
options.IDA.shortCircuit = true;
% SF2 = [0:0.25:1.5 , 2:0.5:5 , 5.75:0.75:8]; % Scale factors to use for each IDA curve
