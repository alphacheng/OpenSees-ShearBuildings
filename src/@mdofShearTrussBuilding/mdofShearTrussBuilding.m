classdef mdofShearTrussBuilding < OpenSeesAnalysis

properties

% General settings

    g                               % Acceleration due to gravity
    units                           % Units used for labels
    verbose = true;                 % Print updates on progress to MATLAB console
    pathOf                          % Paths to OpenSees, external functions, etc.

% Building definition

    nStories                        % Number of stories
    includeExplicitPDelta = true    % Set whether to include gravity loads

    storyMass                       % Mass of each story
    storyHeight                     % Height of each story
    storySpringDefinition           % OpenSees uniaxial material definition for the story springs
    storyTrussDefinition            % OpenSees uniaxial material definition for the trusses

    fundamentalPeriod               % Fundamental period of the structure

    damping_ModeA  = 1;             % Mode A for rayleigh damping
    damping_ModeB  = 3;             % Mode B for rayleigh damping
    damping_RatioA = 0.02;          % Damping ratio for mode A
    damping_RatioB = 0.02;          % Damping ratio for mode B

% Gravity load options

    % optionsGravityLoads - struct containing settings for response history analysis
    %
    % Contains the following fields:
    %   constraints - struct containing settings for the constraints
    %       type - specifies the constraint type to use; select from 'Plain', 'Transformation', and 'Penalty'
    %       penalty - struct containing settings for the penalty method
    %           alphaS - penalty value on single-point constraints
    %           alphaM - penalty value on multi-point constraints
    %
    %   test - struct containing settings for the test method
    %       type        - specifies the test method to use; select from 'NormDispIncr' and 'EnergyIncr'
    %       tolerance   - tolerance for the test
    %       iterations  - max iterations
    %       print       - print flag
    %       normType    - norm type
    %
    %   algorithm - cell vector containing the names of algorithms to be used; select from 'Newton', 'KrylovNewton', and 'ModifiedNewton'
    %
    optionsGravityLoads = struct('constraints' , struct('type'       , 'Transformation', ...
                                                        'penalty'    , struct('alphaS',1.0e12, ...
                                                                              'alphaM',1.0e12) ...
                                 ), ...
                                 'test'        , struct('type'       , 'NormDispIncr', ...
                                                        'tolerance'  , 1e-6, ...
                                                        'iterations' , 10, ...
                                                        'print'      , 1, ...
                                                        'normType'   , 2 ...
                                 ), ...
                                 'algorithm'   , 'Newton' ...
    );

% Pushover analysis options

    % optionsPushover - struct containing settings for pushover analysis
    %
    % Contains the following fields:
    %   controlStory - Control story for the pushover analysis
    %
    %   stepSize - Step size for the pushover analysis
    %
    %   maxDrift - Pushover analysis will abort if the drift of the control story reaches this value
    %
    %   constraints - struct containing settings for the constraints
    %       type - specifies the constraint type to use; select from 'Plain', 'Transformation', and 'Penalty'
    %       penalty - struct containing settings for the penalty method
    %           alphaS - penalty value on single-point constraints
    %           alphaM - penalty value on multi-point constraints
    %
    %   test - struct containing settings for the test method
    %       type        - specifies the test method to use; select from 'NormDispIncr' and 'EnergyIncr'
    %       tolerance   - tolerances for the tests; if switching algorithms fails, analysis will cycle through these
    %       iterations  - max iterations
    %       print       - print flag
    %       normType    - norm type
    %
    %   algorithm - cell vector containing the names of algorithms to be used; select from 'Newton', 'KrylovNewton', and 'ModifiedNewton'
    %
    optionsPushover = struct('controlStory', 'roof', ...
                             'stepSize'    , 0.001, ...
                             'maxDrift'    , 6.0, ...
                             'constraints' , struct('type'       , 'Plain', ...
                                                    'penalty'    , struct('alphaS',1.0e12, ...
                                                                          'alphaM',1.0e12) ...
                            ), ...
                             'test'        , struct('type'       , 'NormDispIncr', ...
                                                    'tolerance'  , [1e-5,1e-4,1e-3], ...
                                                    'iterations' , 10, ...
                                                    'print'      , 1, ...
                                                    'normType'   , 2 ...
                            ), ...
                             'algorithm'   , {{ 'Newton','KrylovNewton','ModifiedNewton' }} ...
    );

% Response history options

    % optionsResponseHistory - struct containing settings for response history analysis
    %
    % Contains the following fields:
    %   constraints - struct containing settings for the constraints
    %       type - specifies the constraint type to use; select from 'Plain', 'Transformation', and 'Penalty'
    %       penalty - struct containing settings for the penalty method
    %           alphaS - penalty value on single-point constraints
    %           alphaM - penalty value on multi-point constraints
    %
    %   test - struct containing settings for the test method
    %       type        - specifies the test method to use; select from 'NormDispIncr' and 'EnergyIncr'
    %       tolerance   - tolerances for the tests; if switching algorithms fails, analysis will cycle through these
    %       iterations  - max iterations
    %       print       - print flag
    %       normType    - norm type
    %
    %   algorithm - cell vector containing the names of algorithms to be used; select from 'Newton', 'KrylovNewton', and 'ModifiedNewton'
    %
    optionsResponseHistory = struct('constraints' , struct('type'       , 'Transformation', ...
                                                           'penalty'    , struct('alphaS',1.0e12, ...
                                                                                 'alphaM',1.0e12) ...
                                    ), ...
                                    'test'        , struct('type'       , 'NormDispIncr', ...
                                                           'tolerance'  , [1e-5,5e-5,1e-4], ...
                                                           'iterations' , 10, ...
                                                           'print'      , 1, ...
                                                           'normType'   , 2 ...
                                    ), ...
                                    'algorithm'   , {{ 'Newton','KrylovNewton','ModifiedNewton' }} ...
    );

% Incremental dynamic analysis options

    SNRT    % Geometric mean of the 5% damped spectral intensity at the fundamental period

    % optionsIDA - struct containing settings for incremental dynamic analysis
    %
    % Contains the following fields:
    %   tExtra              - Extra time to add to end of analysis
    %   nMotions            - Number of ground motions to analyze
    %   ST                  - Vector of intensities to scale each ground motion to
    %   collapseDriftRatio  - Story drift ratio used to define collapse
    %   collapseProbability - Collapse probability being assessed
    %   rating_DR           - Qualitative rating of the design requirements
    %   rating_TD           - Qualitative rating of the test data
    %   rating_MDL          - Qualitative rating of the archetype models
    %
    optionsIDA = struct('tExtra',5, ...
                        'nMotions',7, ...
                        'ST',0.25:0.25:8, ...
                        'collapseDriftRatio',0.05, ...
                        'rating_DR','C', ...
                        'rating_TD','C', ...
                        'rating_MDL','C', ...
                        'shortCircuit',true, ...
                        'ST_tol',0.1, ...
                        'ST_step',0.5 ...
    );

% Equivalent Lateral Force options

    seismicDesignCategory       % Seismic design category (ASCE 7-10 Section 11.6)
    impFactor                   % Seismic importance factor (ASCE 7-10 Section 1.5)
    respModCoeff                % Response modification coefficient (ASCE 7-10 Section 12.2)
    deflAmplFact                % Deflection amplification factor (ASCE 7-10 Section 12.2)
    overstrengthFactor          % Overstrength factor (ASCE 7-10 Section 12.2)

end
properties (Access = protected)
    validAlgorithms  = {'Newton','KrylovNewton','ModifiedNewton'};
    validConstraints = {'Plain','Penalty','Transformation'};
    validTests       = {'NormDispIncr','EnergyIncr'};
end %properties

methods

%###############################################################################
%% Constructor and set methods #################################################
%###############################################################################

function obj = mdofShearTrussBuilding(nStories)
    obj.nStories = nStories;
end

%###############################################################################
%% Shared model functions ######################################################
%###############################################################################

function constructBuilding(obj,fid)
%% CONSTRUCTBUILDING Create the OpenSees model based on current properties
%
%    CONSTRUCTBUILDING(obj,fid) writes the OpenSees code that represents
%       the model to the file specified by fid.
%

    fprintf(fid,'source [file join {%s} {updateRayleighDamping.tcl}]\n',obj.pathOf.tclfunctions);
    fprintf(fid,'model BasicBuilder -ndm 2 -ndf 2 \n');
    % Nodes ----------------------------------------------------------------
    fprintf(fid,'node 0 0.0 0.0\n');
    for i = 1:obj.nStories
        height = sum(obj.storyHeight(1:i));
        fprintf(fid,'node %i 0.0 %g\n',i,height);
        fprintf(fid,'node %i 0.0 %g\n',i*10,height);
    end
    % Masses ---------------------------------------------------------------
    for i = 1:obj.nStories
        fprintf(fid,'mass %i %g %g \n',i*10,obj.storyMass(i),obj.storyMass(i));
    end
    % Fixity ---------------------------------------------------------------
    fprintf(fid,'fix 0 1 1\n');
    for i = 1:obj.nStories
        fprintf(fid,'equalDOF %i %i 1 2\n',(i-1)*10,i);
    end
    % Materials ------------------------------------------------------------
    for i = 1:length(obj.storySpringDefinition)
        fprintf(fid,'%s\n',obj.storySpringDefinition{i});
    end
    for i = 1:length(obj.storyTrussDefinition)
        fprintf(fid,'%s\n',obj.storyTrussDefinition{i});
    end
    % Elements -------------------------------------------------------------
    for i = 1:obj.nStories
        fprintf(fid,'element zeroLength %i %i %i -mat %i -dir 1 -doRayleigh 1\n',i,i,i*10,i);
    end
    for i = 1:obj.nStories
        fprintf(fid,'element corotTruss %i %i %i 1 %i\n',i*10,(i-1)*10,i*10,i*10);
    end

end %function:constructBuilding

function applyGravityLoads(obj,fid)
%% APPLYGRAVITYLOADS Write gravity load commands to file

    fprintf(fid,'pattern Plain 0 Linear {\n');
    for i = 1:obj.nStories
        fprintf(fid,'    load %i 0 -%g\n',i*10,obj.storyMass(i)*obj.g);
    end
    fprintf(fid,'}\n');
    fprintf(fid,'system UmfPack\n');
    switch obj.optionsGravityLoads.constraints.type
        case 'Penalty'
            fprintf(fid,'constraints Penalty %g %g\n',obj.optionsGravityLoads.constraints.penalty.alphaS,...
                                                      obj.optionsGravityLoads.constraints.penalty.alphaM);
        otherwise
            fprintf(fid,'constraints %s\n',obj.optionsGravityLoads.constraints.type);
    end
    fprintf(fid,'numberer RCM\n');
    fprintf(fid,'test %s %g %i %i %i\n',obj.optionsGravityLoads.test.type,...
                                        obj.optionsGravityLoads.test.tolerance,...
                                        obj.optionsGravityLoads.test.iterations,...
                                        obj.optionsGravityLoads.test.print,...
                                        obj.optionsGravityLoads.test.normType);
    fprintf(fid,'algorithm %s\n',obj.optionsGravityLoads.algorithm);
    fprintf(fid,'integrator LoadControl 0.1\n');
    fprintf(fid,'analysis Static\n');
    fprintf(fid,'analyze 10\n');
    fprintf(fid,'loadConst -time 0.0\n');
    fprintf(fid,'wipeAnalysis\n');
end %function:applyGravityLoads


%###############################################################################
%% Analyses ####################################################################
%###############################################################################

function [eigenvals,eigenvecs] = eigenvalues(obj)
%% EIGENVALUES Eigenvalue analysis of system.
%
%   eigenvals = EIGENVALUES(obj) returns the eigenvalues of obj in
%       the vector eigenvals. Note that the eigenvalues are equal to
%       the square of the circular natural frequencies, not the
%       frequencies themselves.
%
%   [eigenvals,eigenvecs] = EIGENVALUES(obj) returns the eigenvalues
%       in the vector eigenvals and the eigenvectors of the first
%       mode of all nodes of obj in the vector eigenvecs.
%

    filename_input = obj.scratchFile('mdofShearBuilding2d_input.tcl');
    filename_vals  = obj.scratchFile('mdofShearBuilding2d_vals.out');
    filename_vecs  = obj.scratchFile('mdofShearBuilding2d_vecs.out');

    fid = fopen(filename_input,'w');

    obj.constructBuilding(fid)
    if obj.includeExplicitPDelta
        obj.applyGravityLoads(fid)
    end

    fprintf(fid,'set eigs [eigen -fullGenLapack %i]\n',obj.nStories);
    fprintf(fid,'set eigenvalues $eigs\n');
    fprintf(fid,'set vecs {}\n');
    fprintf(fid,'set vecfid [open %s w+]\n',obj.path_for_tcl(filename_vecs));
    fprintf(fid,'for {set i 1} {$i <= %i} {incr i} {\n',obj.nStories);
    fprintf(fid,'  lappend vecs [nodeEigenvector [expr $i*10] 1]\n');
    fprintf(fid,'  puts $vecfid [lindex $vecs [expr $i - 1] 0]\n');
    fprintf(fid,'}\n');
    fprintf(fid,'set eigfid [open %s w+]\n',obj.path_for_tcl(filename_vals));
    fprintf(fid,'puts $eigfid $eigs\n');
    fprintf(fid,'close $eigfid\n');
    fprintf(fid,'close $vecfid\n');
    fclose(fid);

    [~,~] = obj.runOpenSees(filename_input);

    eigenvals = dlmread(filename_vals);
    if nargout == 2
        eigenvecs = dlmread(filename_vecs);
    end

    if obj.deleteFilesAfterAnalysis
        delete(filename_input,filename_vals,filename_vecs);
    end

end %function:eigenvalues

function results = responseHistory(obj,gmFile,dt,SF,tEnd,gmID,indexNum)
%% RESPONSEHISTORY Perform response history analysis
%
%   results = RESPONSEHISTORY(obj,gmFile,dt,SF,tEnd,gmID,indexNum)
%       Returns the results of a response history analysis of obj
%       subject to ground motion stored in gmFile with timestep dt
%       scaled by SF. Analysis concludes at tEnd. gmID and indexNum
%       are used for incremental dynamic analyses and are optional
%       if an IDA is not being conducted.
%
%   results has the following fields:
%
%   gmID         -
%   indexNum     -
%   SF           - Scale factor used in the analysis
%   textOutput   - Text output from OpenSees
%   groundMotion - Scaled ground motion used as input
%   displacement   - Time history of the total drift of each story
%   storyShear   - Time history of story shears
%   storyDrift   - Time history of story drifts
%   roofDrift    - Time history of the total roof drift
%   baseShear    - Time history of the base shear
%

% Initialize Results
    results = struct;
    if nargin < 6
        gmID = '01a';
        indexNum = 1;
    else
        results.gmID = gmID;
        results.indexNum = indexNum;
    end
    results.SF = SF;

    % Filenames
    filenames.input              = obj.scratchFile(sprintf('mdofShearTrussBuilding_input_%s_%i.tcl',gmID,indexNum));
    filenames.output_timeSeries  = obj.scratchFile(sprintf('mdofShearTrussBuilding_timeSeries_%s_%i.out',gmID,indexNum));
    filenames.output_def_x       = obj.scratchFile(sprintf('mdofShearTrussBuilding_disp_x_%s_%i.out',gmID,indexNum));
    filenames.output_def_y       = obj.scratchFile(sprintf('mdofShearTrussBuilding_disp_y_%s_%i.out',gmID,indexNum));
    filenames.output_vel_x       = obj.scratchFile(sprintf('mdofShearTrussBuilding_vel_x_%s_%i.out',gmID,indexNum));
    filenames.output_vel_y       = obj.scratchFile(sprintf('mdofShearTrussBuilding_vel_y_%s_%i.out',gmID,indexNum));
    filenames.output_force_story = obj.scratchFile(sprintf('mdofShearTrussBuilding_force_s_%s_%i.out',gmID,indexNum));
    filenames.output_force_truss = obj.scratchFile(sprintf('mdofShearTrussBuilding_force_t_%s_%i.out',gmID,indexNum));

    % Create .tcl file
    fid = fopen(filenames.input,'w');

    obj.constructBuilding(fid)
    if obj.includeExplicitPDelta
        obj.applyGravityLoads(fid)
    end

    fprintf(fid,'timeSeries Path 1 -dt %g -filePath {%s} -factor %g\n',dt,gmFile,SF);
    fprintf(fid,'pattern UniformExcitation 1 1 -accel 1\n');

    fprintf(fid,'recorder Node -file {%s} -time -timeSeries 1 -node 0 -dof 1 accel\n',filenames.output_timeSeries);
    fprintf(fid,'recorder Node -file {%s} -nodeRange 10 %i -dof 1 disp \n',filenames.output_def_x,obj.nStories*10);
    fprintf(fid,'recorder Node -file {%s} -nodeRange 10 %i -dof 2 disp \n',filenames.output_def_y,obj.nStories*10);
    fprintf(fid,'recorder Node -file {%s} -nodeRange 10 %i -dof 1 vel  \n',filenames.output_vel_x,obj.nStories*10);
    fprintf(fid,'recorder Node -file {%s} -nodeRange 10 %i -dof 2 vel  \n',filenames.output_vel_y,obj.nStories*10);
    fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i force \n',filenames.output_force_story,obj.nStories);
    fprintf(fid,'recorder Element -file {%s} -eleRange 10 %i force \n',filenames.output_force_truss,obj.nStories*10);
    fprintf(fid,'record \n');

    fprintf(fid,'system UmfPack \n');
    switch obj.optionsResponseHistory.constraints.type
        case 'Penalty'
            fprintf(fid,'constraints Penalty %g %g\n',obj.optionsResponseHistory.constraints.penalty.alphaS,...
                                                      obj.optionsResponseHistory.constraints.penalty.alphaM);
        otherwise
            fprintf(fid,'constraints %s\n',obj.optionsResponseHistory.constraints.type);
    end
    % test
    testArgs = cell(length(obj.optionsResponseHistory.test.tolerance),1);
    for i = 1:length(obj.optionsResponseHistory.test.tolerance)
        testArgs{i} = sprintf('%s %g %i %i %i',obj.optionsResponseHistory.test.type,...
                                               obj.optionsResponseHistory.test.tolerance(i),...
                                               obj.optionsResponseHistory.test.iterations,...
                                               obj.optionsResponseHistory.test.print,...
                                               obj.optionsResponseHistory.test.normType);
    end

    fprintf(fid,'numberer RCM \n');

    fprintf(fid,'updateRayleighDamping %i %g %i %g\n',...
        obj.damping_ModeA,obj.damping_RatioA,...
        obj.damping_ModeB,obj.damping_RatioB);

    fprintf(fid,'integrator Newmark 0.50 0.25\n');
    fprintf(fid,'analysis VariableTransient \n');

    fprintf(fid,'set currentTime [getTime]\n');
    fprintf(fid,'while { $currentTime < %g } {\n',tEnd);
    fprintf(fid,'    algorithm %s\n',obj.optionsResponseHistory.algorithm{1});
    fprintf(fid,'    test %s\n',testArgs{1});
    fprintf(fid,'    set ok [analyze 1 %g]\n',dt);
    for i = 1:length(testArgs)
        if i == 1; k = 2; else; k = 1; end
        for j = k:length(obj.optionsResponseHistory.algorithm)
            fprintf(fid,'    if { $ok != 0 } {\n');
            fprintf(fid,'        algorithm %s\n',obj.optionsResponseHistory.algorithm{j});
            fprintf(fid,'        test %s\n',testArgs{i});
            fprintf(fid,'        set ok [analyze 1 %g]\n',dt);
            fprintf(fid,'    }\n');
        end
    end
    fprintf(fid,'    if { $ok != 0 } {\n');
    fprintf(fid,'        exit 2\n');
    fprintf(fid,'    }\n');
    fprintf(fid,'    set currentTime [getTime]\n');
    fprintf(fid,'}\n');

    fprintf(fid,'exit 1 \n');
    fclose(fid);

    % Run OpenSees
    [status, result] = obj.runOpenSees(filenames.input);
    results.textOutput = result;
    switch status
        case 1
            results.exitStatus = 'Analysis Successful';
        case 2
            results.exitStatus = 'Analysis Failed';
        otherwise
            fprintf('%s\n',result);
            error('Analysis Failed in Unknown Manner (exit code: %i)',status);
    end

    % Read Results
    temp = dlmread(filenames.output_timeSeries);
    results.time = temp(:,1);
    results.groundMotion = temp(:,2);
    results.displacement_x = dlmread(filenames.output_def_x);
    results.displacement_y = dlmread(filenames.output_def_y);
    results.velocity_x = dlmread(filenames.output_vel_x);
    results.velocity_y = dlmread(filenames.output_vel_y);
    temp = dlmread(filenames.output_force_story);
    results.storyShear = temp(:,1:4:end);
    temp = dlmread(filenames.output_force_truss);
    results.trussForce = zeros(size(temp,1),size(temp,2)/2);
    results.trussForce(:,1:2:end) = temp(:,1:4:end);
    results.trussForce(:,2:2:end) = temp(:,2:4:end);

    % Computed Results
    storyDrift = results.displacement_x;
    storyDrift(:,2:end) = storyDrift(:,2:end)-storyDrift(:,1:(end-1));
    results.storyDrift = storyDrift;
    results.roofDrift = results.displacement_x(:,end);
    results.baseShear = results.storyShear(:,1);

    % Clean Folder
    if obj.deleteFilesAfterAnalysis
        fields = fieldnames(filenames);
        for i = 1:length(fields)
            delete(filenames.(fields{i}));
        end
    end
end %function:responseHistory

function energy = energyCriterion(obj,results)
%% ENERGYCRITERION Calculate the energy collapse criterion

    time = results.time;

    M = obj.storyMass(:);

    u_dot     = results.velocity_x;
    u_ddot_eq = results.groundMotion;

    expression = u_dot*M.*u_ddot_eq;

    E_EQ = -cumtrapz(time,expression);

    u = results.displacement_y;
    E_G = -u*M*obj.g;
    E_G_norm = E_G - E_G(1);

    collapseIndex = find(E_G_norm > E_EQ, 1);
    if isempty(collapseIndex)
        collapseIndex = NaN;
        collapse = false;
    else
        collapse = true;
    end

    energy = struct;
    energy.collapse = collapse;
    energy.collapseIndex = collapseIndex;
    energy.earthquake   = E_EQ;
    energy.gravity      = E_G;
    energy.norm_gravity = E_G_norm;
end

function results = incrementalDynamicAnalysis(obj,gm_mat)
%% INCREMENTALDYNAMICANALYSIS Perform an IDA on the model
%
%   results = incrementalDynamicAnalysis(obj,gm_mat) uses the .mat file specified
%       by `gm_mat` to perform an incremental dynamic analysis. The .mat file must
%       be of the format generated by the `FEMAP695` package.
%

    SMT = FEMAP695_SMT(obj.fundamentalPeriod,obj.seismicDesignCategory);
    SF1 = FEMAP695_SF1(obj.fundamentalPeriod,obj.seismicDesignCategory);

    temp = load(gm_mat);
    ground_motions = temp.ground_motions; % necessary for parfor transparency

    gm = cell(obj.optionsIDA.nMotions,1);
    gm(:) = {struct};
    parfor gmIndex = 1:obj.optionsIDA.nMotions
        totaltic = tic;
        if obj.verbose %#ok obviously the class object needs to be broadcast? Yipes
            w = getCurrentWorker();
            fprintf('Calculating ground motion %i on %i\n',gmIndex,w.ProcessId)
        end
        gmID = ground_motions(gmIndex).ID;
        gmFile = scratchFile(obj,sprintf('acc%s.acc',gmID));
        dlmwrite(gmFile,ground_motions(gmIndex).normalized_acceleration*obj.g);
        dt     = ground_motions(gmIndex).dt;
        tEnd   = ground_motions(gmIndex).time(end) + obj.optionsIDA.tExtra;

        % solve_tic = tic;
        % fun = @(SF) IDApoint(obj,gmFile,dt,SF,tEnd,gmID,1) - 1;
        % fzero_options = optimset(optimset('fzero'),'TolX',obj.optionsIDA.ST_tol);
        % collapseSF = fzero(fun,3,fzero_options);
        % while fun(collapseSF) < 0
        %     collapseSF = collapseSF + obj.optionsIDA.ST_tol;
        % end
        % solvetime = toc(solve_tic);
        % collapseST = collapseSF*SMT/SF1;
        % ST = [obj.optionsIDA.ST_step:obj.optionsIDA.ST_step:(collapseST-2*obj.optionsIDA.ST_step) (collapseST-obj.optionsIDA.ST_step):obj.optionsIDA.ST_step/10:collapseST];
        % SF = ST*SF1/SMT;

        % calc_tic = tic;
        % nHistories = length(SF);
        % rh = cell(nHistories,1);
        % for index = 1:nHistories
        %     rh{index} = obj.responseHistory(gmFile,dt,SF(index),tEnd,gmID,index);
        %     rh{index}.energy = obj.energyCriterion(rh{index});
        %     rh{index}.maxDriftRatio = max(max(abs(rh{index}.storyDrift))./obj.storyHeight);
        % end
        % rh = [rh{:}];
        % calctime = toc(calc_tic);

        % fun = @(ST) IDApoint(obj,gmFile,dt,ST*SF1/SMT,tEnd,gmID) - 1;
        ST = 0;
        count = 0;
        ER = 0;
        ST_all = NaN(ceil(3/obj.optionsIDA.ST_step)*2,1);
        rh = cell(size(ST_all));
        while ER < 1
            ST = ST + obj.optionsIDA.ST_step;
            count = count + 1;
            ST_all(count) = ST;
            [ER,rh{count}] = IDApoint(obj,gmFile,dt,ST*SF1/SMT,tEnd,gmID,1);
        end
        ST_collapse = ST;
        ST_noCollapse = ST_all(end-1);
        if isnan(ST_noCollapse)
            ST_noCollapse = 0;
        end
        % Collapse! Now to bisect
        while (ST_collapse - ST_noCollapse) > obj.optionsIDA.ST_tol
            if ER > 1 %Collapse occurred, need to drop
                ST_collapse = ST;
            else %Collapse didn't occur, need to step up
                ST_noCollapse = ST;
            end
            ST = ST_noCollapse + (ST_collapse - ST_noCollapse)/2;
            count = count + 1;
            ST_all(count) = ST;
            [ER,rh{count}] = IDApoint(obj,gmFile,dt,ST*SF1/SMT,tEnd,gmID,1);
        end
        [ST_all,I] = sort(ST_all);
        ST_all(isnan(ST_all)) = [];
        rh = rh(I);
        rh = [rh{:}];

        temp = [rh.energy];
        gm{gmIndex}.ID = gmID;
        gm{gmIndex}.maxDriftRatio = [rh.maxDriftRatio];
        gm{gmIndex}.collapseIndex = find(~isnan([temp.collapseIndex]),1);
        gm{gmIndex}.ST = ST_all';
        gm{gmIndex}.SCT = ST_all(gm{gmIndex}.collapseIndex);
        gm{gmIndex}.rh = rh;
        % gm{gmIndex}.solvetime = solvetime;
        % gm{gmIndex}.calctime  = calctime;

        if obj.deleteFilesAfterAnalysis
            delete(gmFile)
        end
        if obj.verbose
            fprintf('Finished calculating ground motion %i on process %i in %g seconds\n',gmIndex,w.ProcessId,toc(totaltic))
        end
    end
    gm = [gm{:}];
    SCT_hat = median([gm.SCT]);

    results = struct;
    results.gm      = gm;
    results.SMT     = SMT;
    results.SCT_hat = SCT_hat;

end

function [energyRatio,results] = IDApoint(obj,gmFile,dt,SF,tEnd,gmID,index)
% IDAPOINT Perform a response history and return the maximum gravity to earthquake energy ratio.
%
%   This function is only used by INCREMENTALDYNAMICANALYSIS, but must be
%       separate for parallelization reasons.
%
    results = responseHistory(obj,gmFile,dt,SF,tEnd,gmID,index);
    results.energy = energyCriterion(obj,results);
    results.maxDriftRatio = max(max(abs(results.storyDrift))./obj.storyHeight);
    energyRatio = max(results.energy.norm_gravity ./ (results.energy.earthquake + eps)); % add eps to earthquake to avoid 0/0

end

%###############################################################################
%% Design Stuff ################################################################
%###############################################################################

function results = ELFanalysis(obj)
%% ELFANALYSIS Equivalent Lateral Force procedure (ASCE 7-10)
%
%   results = ELFANALYSIS(obj) generates the design story forces,
%       shears, and allowable story drifts for the information in obj.
%
%   The struct results contains the following fields:
%
%   seismicResponseCoefficient      ; C_s (Section 12.8.1.1)
%   baseShear                       ; V (Section 12.8.1)
%   storyForce                      ; F_x (Section 12.8.3)
%   storyShear                      ; V_x (Section 12.8.4)
%   allowableDrift                  ; Delta_a (Table 12.12-1)
%

    results = struct;

    SDS = FEMAP695_mappedValue('SDS',obj.seismicDesignCategory);
    SD1 = FEMAP695_mappedValue('SD1',obj.seismicDesignCategory);

    approxFundamentalPeriod = 0.02*sum(obj.storyHeight)^0.75;

    if SD1 <= 0.1
        Cu = 1.7;
    elseif SD1 >= 0.4
        Cu = 1.4;
    else
        Cu = interp1([0.1 0.15 0.2 0.3 0.4],[1.7 1.6 1.5 1.4 1.4],SD1);
    end

    if isempty(obj.fundamentalPeriod)
        obj.fundamentalPeriod = 0.02*sum(obj.storyHeight)^0.75;
    elseif obj.fundamentalPeriod > Cu*approxFundamentalPeriod
        obj.fundamentalPeriod = Cu*approxFundamentalPeriod;
    end

    maxSeismicResponseCoefficient = SD1/(obj.fundamentalPeriod*obj.respModCoeff/obj.impFactor);
    results.seismicResponseCoefficient = min(SDS/(obj.respModCoeff/obj.impFactor),maxSeismicResponseCoefficient);

    seismicWeight = sum(obj.storyMass)*obj.g;
    results.baseShear = seismicWeight*results.seismicResponseCoefficient;

    if obj.fundamentalPeriod <= 0.5
        k = 1;
    elseif obj.fundamentalPeriod >= 2.5
        k = 2;
    else
        k = interp1([0.5 2.5],[1 2],obj.fundamentalPeriod);
    end

    verticalDistributionFactor = (obj.storyMass*obj.g .* cumsum(obj.storyHeight).^k)/sum(obj.storyMass*obj.g .* cumsum(obj.storyHeight).^k);

    results.storyForce = verticalDistributionFactor*results.baseShear;
    results.storyShear = zeros(1,obj.nStories);
    for i = 1:obj.nStories
        results.storyShear(i) = sum(results.storyForce(i:end));
    end

    results.allowableDrift = 0.020*obj.storyHeight;

end %function:ELFanalysis

function spring = springDesign(obj,analysisResults,springGivens)
%% SPRINGDESIGN  Story spring design

    spring = struct;

    designStiffness = springGivens.stiffnessSafety*obj.deflAmplFact*analysisResults.storyShear./(obj.impFactor*analysisResults.allowableDrift);

    designStrength  = springGivens.strengthSafety*obj.overstrengthFactor*analysisResults.storyShear;

    if springGivens.enforceMinimumStiffness
        for i = 2:length(designStiffness)
            if designStiffness(i) < springGivens.minimumRatio*designStiffness(i-1)
                designStiffness(i) = springGivens.minimumRatio*designStiffness(i-1);
            end
        end
    end
    if springGivens.enforceMinimumStrength
        for i = 2:length(designStrength)
            if designStrength(i) < springGivens.minimumRatio*designStrength(i-1)
                designStrength(i) = springGivens.minimumRatio*designStrength(i-1);
            end
        end
    end

    for i = 1:obj.nStories
        if springGivens.includePDelta
            Px = sum(obj.storyMass(i:end))*obj.g;
            theta = (Px*analysisResults.allowableDrift(i)*obj.impFactor)/(analysisResults.storyShear(i)*obj.storyHeight(i)*obj.deflAmplFact);
            theta_max = min(0.5/obj.deflAmplFact,0.25);
            theta = min(theta,theta_max);
        else
            theta = 0;
        end
        spring(i).K0       = (1-theta)*designStiffness(i);  % elastic stiffness
        spring(i).as       = springGivens.as - theta;       % strain hardening ratio
        spring(i).ad       = springGivens.ad + theta;       % strain hardening ratio
        spring(i).Lambda_S = springGivens.Lambda_S;         % Cyclic deterioration parameter - strength
        spring(i).Lambda_K = springGivens.Lambda_K;         % Cyclic deterioration parameter - stiffness
        spring(i).c_S      = springGivens.c_S;              % rate of deterioration - strength
        spring(i).c_K      = springGivens.c_K;              % rate of deterioration - stiffness
        spring(i).Res      = springGivens.Res;              % residual strength ratio
        spring(i).D        = springGivens.D;                % rate of cyclic deterioration
        spring(i).nFactor  = springGivens.nFactor;          % elastic stiffness amplification factor
        spring(i).theta    = theta;

        spring(i).V_c = (1-theta)*designStrength(i);        % strength at capping
        spring(i).V_y = springGivens.C_yc*spring(i).V_c;    % effective yield strength

        spring(i).defl_y  = spring(i).V_y./spring(i).K0;                                                    % deflection at yield
        spring(i).defl_p  = (spring(i).V_c-spring(i).V_y)./((spring(i).as+theta)*spring(i).K0);             % pre-capping deflection
        spring(i).defl_pc = spring(i).V_c/(spring(i).ad*spring(i).K0);                                      % post-capping deflection
        spring(i).defl_u  = springGivens.C_upc*(spring(i).defl_y + spring(i).defl_p + spring(i).defl_pc);   % ultimate deflection capacity

        spring(i).definition = bilinearMaterialDefinition(i,spring(i));
    end
end %function:springDesign

%###############################################################################
%% Plot Functions ##############################################################
%###############################################################################

function animateResponseHistory(obj,results,dt)

    if nargin == 2
        dt = max(diff(results.time));
    end
    cumHeights = cumsum(obj.storyHeight);

    figure

    grid on
    grid minor

    xMax = max(max(abs(results.displacement_x)));
    yMax = sum(obj.storyHeight) + obj.storyHeight(1);

    axis([-xMax xMax 0 yMax]);
    ax = gca;
    ax.YTick = cumHeights;

    ylabel(sprintf('Height (%s)',obj.units.length))
    xlabel(sprintf('Drift (%s)',obj.units.length))

    xPos = [0 results.displacement_x(1,:)];
    yPos = [0 cumHeights] + [0 results.displacement_y(1,:)];

    h = animatedline(xPos,yPos,'Marker','*');
    t = text(2/3*xMax,1/5*yMax,'Time: 0.0s');

    for i = 1:size(results.displacement_x,1)
        displayText = sprintf('Time: %4.1fs',i*dt);
        delete(t)
        clearpoints(h)
        xPos = [0 results.displacement_x(i,:)];
        yPos = [0 cumHeights] + [0 results.displacement_y(i,:)];
        addpoints(h,xPos,yPos)
        t = text(ax,2/3*xMax,1/5*yMax,displayText);
        drawnow
    end
end


end %methods

end %classdef
