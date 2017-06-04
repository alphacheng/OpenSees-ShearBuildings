classdef mdofShearBuilding2d < OpenSeesAnalysis
% MDOFSHEARBUILDING2D

    properties
    % General settings

        g                               % Acceleration due to gravity
        units                           % Units used for labels

    % Building definition

        nStories                        % Number of stories

        storyMass                       % Mass of each story
        storyHeight                     % Height of each story
        storyStiffness                  % Elastic stiffness of each story
        storySpringDefinition           % OpenSees uniaxial material definition

        fundamentalPeriod               % Fundamental period of the structure

    % Pushover analysis options

        controlStory = 'roof';          % Control story for the pushover analysis
        pushover_stepSize   = 0.001;    % Step size for the pushover analysis
        pushover_maxDrift   = 6.0;      % Pushover analysis will abort if the drift of the control story reaches this value

        % optionsPushover - struct containing settings for pushover analysis
        %   Defaults for optionsPushover are set by the constructor.
        %
        % Contains the following fields:
        %   constraints - struct containing settings for the constraints
        %       type - specifies the constraint type to use; select from 'Plain' and 'Penalty'
        %       penalty - struct containing settings for the penalty method
        %           alphaS -
        %           alphaM -
        %
        %   test - struct containing settings for the test method
        %       type        - specifies the test method to use; select from 'NormDispIncr' and 'EnergyIncr'
        %       tolerance   - tolerances for the tests; if switching algorithms fails, analysis will cycle through these
        %       iterations  - max iterations
        %       print       - print flag
        %       normType    - norm type
        %
        %   algorithm - cell vector containing the names of algorithms to be used
        %
        optionsPushover

    % Response history options

        damping_ModeA  = 1;             % Mode A for rayleigh damping
        damping_ModeB  = 3;             % Mode B for rayleigh damping
        damping_RatioA = 0.02;          % Damping ratio for mode A
        damping_RatioB = 0.02;          % Damping ratio for mode B

        % optionsResponseHistory
        optionsResponseHistory

    % Incremental dynamic analysis options

        % optionsIDA - struct containing settings for incremental dynamic analysis
        %   Defaults for optionsIDA are set by the constructor.
        %
        % Contains the following fields:
        %
        %
        optionsIDA

    % Equivalent Lateral Force options

        seismicDesignCategory       % Seismic design category (Section )
        impFactor                   % Seismic importance factor (Section )
        respModCoeff                % Response modification coefficient (Section )
        deflAmplFact                % Deflection amplification factor (Section )
        overstrengthFactor          % Overstrength factor (Section )

    end
    properties (Access = protected)
        validAlgorithms  = {'Newton','KrylovNewton','ModifiedNewton'};
        validConstraints = {'Plain','Penalty','Transformation'};
        validTests       = {'NormDispIncr','EnergyIncr'};
    end

    methods
        %% Constructor
        function obj = mdofShearBuilding2d(nStories)
            % Constructor
            %
            % obj = mdofShearBuilding2d(nStories)
            %
            obj.nStories = nStories;

            % Pushover options
            obj.optionsPushover.constraints.type = 'Plain';
            obj.optionsPushover.constraints.penalty.alphaS = 1.0e12;
            obj.optionsPushover.constraints.penalty.alphaM = 1.0e12;

            obj.optionsPushover.test.type       = 'NormDispIncr';
            obj.optionsPushover.test.tolerance  = [1e-5,1e-4,1e-3];
            obj.optionsPushover.test.iterations = 10;
            obj.optionsPushover.test.print      = 1;
            obj.optionsPushover.test.normType   = 2;

            obj.optionsPushover.algorithm = { 'Newton','KrylovNewton','ModifiedNewton' };

            % Response history options
            obj.optionsResponseHistory.constraints.type = 'Transformation';
            obj.optionsResponseHistory.constraints.penalty.alphaS = 1.0e12;
            obj.optionsResponseHistory.constraints.penalty.alphaM = 1.0e12;

            obj.optionsResponseHistory.test.type       = 'NormDispIncr';
            obj.optionsResponseHistory.test.tolerance  = [1e-5,5e-5,1e-4];
            obj.optionsResponseHistory.test.iterations = 10;
            obj.optionsResponseHistory.test.print      = 1;
            obj.optionsResponseHistory.test.normType   = 2;

            obj.optionsResponseHistory.algorithm = { 'Newton','KrylovNewton','ModifiedNewton' };

            % Incremental dynamic analysis options
            obj.optionsIDA.tExtra = 5;
            obj.optionsIDA.nMotions = 7;

        end

        %% Set functions
        function set.storyMass(obj,storyMass)
            if ~isnumeric(storyMass)
                error('storyMass should be numeric');
            end
            if isvectorsize(storyMass,obj.nStories) %#ok - nStories is set by the constructor
                obj.storyMass = storyMass;
            else
                error('storyMass should be vector of length %i (number of stories)',obj.nStories); %#ok - nStories is set by the constructor
            end
        end
        function set.storyStiffness(obj,storyStiffness)
            if ~isnumeric(storyStiffness)
                error('storyStiffness should be numeric');
            end
            if isvectorsize(storyStiffness,obj.nStories) %#ok - nStories is set by the constructor
                obj.storyStiffness = storyStiffness;
            else
                error('storyStiffness should be vector of length %i (number of stories)',obj.nStories); %#ok - nStories is set by the constructor
            end
        end
        function set.storySpringDefinition(obj,storySpringDefinition)
            if ~iscell(storySpringDefinition)
                error('storySpringDefinition should be a cell vector');
            end
            if isvectorsize(storySpringDefinition,obj.nStories) %#ok - nStories is set by the constructor
                obj.storySpringDefinition = storySpringDefinition;
            else
                error('storySpringDefinition should be cell vector of length %i (number of stories)',obj.nStories); %#ok - nStories is set by the constructor
            end
        end
        function set.optionsPushover(obj,optionsPushover)
            if isfield(optionsPushover,'constraints')
                optionsPushover.constraints = obj.checkConstraints(optionsPushover.constraints);
            end
            if isfield(optionsPushover,'test')
                optionsPushover.test = obj.checkTest(optionsPushover.test);
            end
            if isfield(optionsPushover,'algorithm')
                optionsPushover.algorithm = obj.checkAlgorithm(optionsPushover.algorithm);
            end
            obj.optionsPushover = optionsPushover;
        end
        function constraints = checkConstraints(obj,constraints)
            if isfield(constraints,'type')
            assert(ischar(constraints.type),'Constraints type must be a character vector');
            check = strcmpi(constraints.type,obj.validConstraints);
            assert(any(check),'Unknown constraints type: %s',constraints.type);
            constraints.type = obj.validConstraints{check};  % Ensure capitalization is correct
            end
        end
        function test = checkTest(obj,test)
            assert(ischar(test.type),'Test method must be a character vector');
            check = strcmpi(test.type,obj.validTests);
            assert(any(check),'Unknown test method: %s',test.type);
            test.type = obj.validTests{check};  % Ensure capitalization is correct
        end
        function algorithm = checkAlgorithm(obj,algorithm)
            assert(iscell(algorithm),'Algorithm list must be a cell vector')
            for i = 1:length(algorithm)
                check = strcmpi(algorithm{i},obj.validAlgorithms);
                assert(any(check),'Unknown test method: %s',algorithm{i});
                algorithm{i} = obj.validAlgorithms{check};  % Ensure capitalization is correct
            end
        end

        %% Functions
        function s = controlStory1(obj)
            if isequal(obj.controlStory,'roof')
                s = obj.nStories;
            else
                s = obj.controlStory;
            end
        end

        function constructBuilding(obj,fid)
            % CONSTRUCTBUILDING Create the OpenSees model based on current properties
            %
            %    CONSTRUCTBUILDING(obj,fid) writes the OpenSees code that represents
            %       the model to the file specified by fid.
            %

            fprintf(fid,'model BasicBuilder -ndm 1 -ndf 1 \n');
            for i = 0:obj.nStories
                fprintf(fid,'node %i 0.0\n',i);
            end
            for i = 1:obj.nStories
                fprintf(fid,'mass %i %g\n',i,obj.storyMass(i));
            end
            fprintf(fid,'fix 0 1 \n');
            for i = 1:length(obj.storySpringDefinition)
                fprintf(fid,'%s\n',obj.storySpringDefinition{i});
            end
            for i = 1:obj.nStories
                fprintf(fid,'element zeroLength %i %i %i -mat %i -dir 1\n',i,i-1,i,i);
            end

        end %function:constructBuilding

        %% Analyses
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

            fprintf(fid,'set eigs [eigen -fullGenLapack %i]\n',obj.nStories);
            fprintf(fid,'set eigenvalues $eigs\n');
            fprintf(fid,'set vecs {}\n');
            fprintf(fid,'set vecfid [open %s w+]\n',obj.path_for_tcl(filename_vecs));
            fprintf(fid,'for {set i 1} {$i <= %i} {incr i} {\n',obj.nStories);
            fprintf(fid,'  lappend vecs [nodeEigenvector $i 1]\n');
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

        function results = pushover(obj,F,type,varargin)
            %% PUSHOVER Perform a pushover analysis
            %
            %   results = PUSHOVER(obj,F,type,typeArg) performs a pushover
            %       analysis with load distribution specified by F and end
            %       condition defined by type and typeArg.
            %
            %   type accepts the following options:
            %       targetDrift
            %       targetPostPeakRatio
            %
            %   results has the following fields:
            %
            %   F                       Force ratios
            %   targetDrift             Target drift for analysis
            %   targetPostPeakRatio     Target post peak ratio
            %   textOutput              Console output from OpenSees
            %   exitStatus              Reports whether analysis was successful
            %   totalDrift              Time history of total drift of stories
            %   storyShear              Time history of story shears
            %   storyDrift              Time history of story drifts
            %   appliedStoryForce       Time history of applied forces
            %   roofDrift               Time history of total roof drift
            %   baseShear               Time history of base shear
            %

            assert(isnumeric(F) & isvectorsize(F,obj.nStories),...
                'F should be a numeric verctor of length %i (number of stories)',obj.nStories);
            if iscolumn(F)
                F = F';
            end

            % Initialize Results
            results = struct;
            results.F = F;

            switch lower(type)
                case 'targetdrift'
                    targetDrift      = varargin{1};

                    results.targetDrift = targetDrift;
                case 'targetpostpeakratio'
                    targetPostPeakRatio = varargin{1};

                    results.targetPostPeakRatio = targetPostPeakRatio;
                otherwise
                    error('Unknown analysis type: %s',type);
            end

            % Filenames
            filename_input          = obj.scratchFile('mdofShearBuilding2d_input.tcl');
            filename_output_def     = obj.scratchFile('mdofShearBuilding2d_disp.out');
            filename_output_force   = obj.scratchFile('mdofShearBuilding2d_force.out');

            % Create .tcl file
            fid = fopen(filename_input,'w');

            obj.constructBuilding(fid)

            fprintf(fid,'timeSeries Linear 1\n');
            fprintf(fid,'pattern Plain 1 1 {\n');
            for i = 1:obj.nStories
                if F(i) ~= 0
                    fprintf(fid,'    load %i %g \n',i,F(i));
                end
            end
            fprintf(fid,'} \n');

            fprintf(fid,'recorder Node -file {%s} -time -nodeRange 1 %i -dof 1 disp \n',filename_output_def,obj.nStories);
            fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i force \n',filename_output_force,obj.nStories);
            fprintf(fid,'record \n');
            fprintf(fid,'system UmfPack \n');
            % constraints
            switch obj.optionsPushover.constraints.type
            case 'Penalty'
                fprintf(fid,'constraints Penalty %g %g\n',obj.optionsPushover.constraints.penalty.alphaS,...
                                                          obj.optionsPushover.constraints.penalty.alphaM);
            otherwise
                fprintf(fid,'constraints %s\n',obj.optionsPushover.constraints.type);
            end
            % test
            testArgs = cell(length(obj.optionsPushover.test.tolerance),1);
            for i = 1:length(obj.optionsPushover.test.tolerance)
                testArgs{i} = sprintf('%s %g %i %i %i',obj.optionsPushover.test.type,...
                                                       obj.optionsPushover.test.tolerance(i),...
                                                       obj.optionsPushover.test.iterations,...
                                                       obj.optionsPushover.test.print,...
                                                       obj.optionsPushover.test.normType);
            end

            fprintf(fid,'numberer RCM \n');
            fprintf(fid,'integrator DisplacementControl %i 1 %g\n',obj.controlStory1,obj.pushover_stepSize);
            fprintf(fid,'analysis Static \n');

            switch lower(type)
                case 'targetdrift'
                    fprintf(fid,'algorithm Newton\n');
                    fprintf(fid,'set ok [analyze %i]\n',ceil(targetDrift/obj.pushover_stepSize));
                    fprintf(fid,'if { $ok != 0 } {\n');
                    fprintf(fid,'    exit 1\n');
                    fprintf(fid,'}\n');
                case 'targetpostpeakratio'
                    fprintf(fid,'set currentLoad [getTime]\n');
                    fprintf(fid,'set maxLoad $currentLoad\n');
                    fprintf(fid,'while { $currentLoad >= [expr %g*$maxLoad] } {\n',targetPostPeakRatio);
                    fprintf(fid,'    algorithm %s\n',obj.optionsPushover.algorithm{1});
                    fprintf(fid,'    test %s\n',testArgs{1});
                    fprintf(fid,'    set ok [analyze 1]\n');
                    for i = 1:length(testArgs)
                        if i == 1; k = 2; else; k = 1; end
                        for j = k:length(obj.optionsPushover.algorithm)
                            fprintf(fid,'    if { $ok != 0 } {\n');
                            fprintf(fid,'        algorithm %s\n',obj.optionsPushover.algorithm{j});
                            fprintf(fid,'        test %s\n',testArgs{i});
                            fprintf(fid,'        set ok [analyze 1]\n');
                            fprintf(fid,'    }\n');
                        end
                    end
                    fprintf(fid,'    if { $ok != 0 } {\n');
                    fprintf(fid,'        exit 2\n');
                    fprintf(fid,'    }\n');
                    fprintf(fid,'    set currentLoad [getTime]\n');
                    fprintf(fid,'    if { $currentLoad > $maxLoad } {\n');
                    fprintf(fid,'        set maxLoad $currentLoad\n');
                    fprintf(fid,'    }\n');
                    fprintf(fid,'    if { [nodeDisp %i 1] > %g } {\n',obj.controlStory1,obj.pushover_maxDrift);
                    fprintf(fid,'        exit 3\n');
                    fprintf(fid,'    }\n');
                    fprintf(fid,'}\n');
                otherwise
                    error('Unknown analysis type: %s',type);
            end

            fprintf(fid,'exit 1 \n');
            fclose(fid);

            % Run OpenSees
            [status, result] = obj.runOpenSees(filename_input);
            results.textOutput = result;
            switch status
                case 1
                    results.exitStatus = 'Analysis Successful';
                case 2
                    results.exitStatus = 'Analysis Failed';
                case 3
                    results.exitStatus = 'Peak Drift Reached';
                otherwise
                    fprintf('%s\n',result);
                    error('Analysis Failed in Unknown Manner (exit code: %i)',status);
            end

            % Read Results
            temp = dlmread(filename_output_def);
            time = temp(:,1);
            results.totalDrift = temp(:,2:end);
            temp = dlmread(filename_output_force);
            results.storyShear = temp(:,2:2:end);

            % Computed Results
            storyDrift = results.totalDrift;
            storyDrift(:,2:end) = storyDrift(:,2:end)-storyDrift(:,1:(end-1));
            results.storyDrift = storyDrift;
            results.appliedStoryForce = time*F;
            results.roofDrift = results.totalDrift(:,end);
            results.baseShear = results.storyShear(:,1);

            % Analysis
            if strcmp(results.exitStatus,'Analysis Successful')
                results.peakShear = max(results.baseShear);
                peakIndex = results.baseShear == results.peakShear;
                results.peakTotalDrift = results.totalDrift(peakIndex,:);
                peakStoryDrift = results.storyDrift(peakIndex,:);

                postPeakIndex = results.roofDrift > results.roofDrift(peakIndex);
                postPeakShear = results.baseShear(postPeakIndex);
                postPeakTotalDrift = results.totalDrift(postPeakIndex,:);
                postPeakStoryDrift = results.storyDrift(postPeakIndex,:);

                results.peak80Shear = 0.8*results.peakShear;
                results.peak80TotalDrift = interp1(postPeakShear,postPeakTotalDrift,results.peak80Shear);
                peak80StoryDrift = interp1(postPeakShear,postPeakStoryDrift,results.peak80Shear);

                results.peakStoryDriftRatio   = peakStoryDrift./obj.storyHeight;
                results.peak80StoryDriftRatio = peak80StoryDrift./obj.storyHeight;
            end

            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_output_def,filename_output_force);
            end
        end %function:pushover

        function results = responseHistory(obj,groundMotionFilename,dt,SF,tend,gmID,indexNum)
            %% RESPONSEHISTORY Perform response history analysis
            %
            %   results = RESPONSEHISTORY(obj,gmFile,dt,SF,tend,gmID,indexNum)
            %       Returns the results of a response history analysis of obj
            %       subject to ground motion stored in gmFile with timestep dt
            %       scaled by SF. Analysis concludes at tend. gmID and indexNum
            %       are used for incremental dynamic analyses.
            %
            %   results has the following fields:
            %
            %   gmID
            %   indexNum
            %   SF
            %   textOutput
            %   groundMotion
            %   totalDrift
            %   storyShear
            %   storyDrift
            %   roofDrift
            %   baseShear
            %

            % Initialize Results
            results = struct;
            results.gmID = gmID;
            results.indexNum = indexNum;
            results.SF = SF;

            % Filenames
            filename_input              = obj.scratchFile(sprintf('mdofShearBuilding2d_input_%s_%i.tcl',gmID,indexNum));
            filename_output_timeSeries  = obj.scratchFile(sprintf('mdofShearBuilding2d_timeSeries_%s_%i.out',gmID,indexNum));
            filename_output_def         = obj.scratchFile(sprintf('mdofShearBuilding2d_disp_%s_%i.out',gmID,indexNum));
            filename_output_force       = obj.scratchFile(sprintf('mdofShearBuilding2d_force_%s_%i.out',gmID,indexNum));

            % Create .tcl file
            fid = fopen(filename_input,'w');

            writeFunction_updateRayleighDamping(fid)
            obj.constructBuilding(fid)

            fprintf(fid,'timeSeries Path 1 -dt %g -filePath {%s} -factor %g\n',dt,groundMotionFilename,SF);
            fprintf(fid,'pattern UniformExcitation 1 1 -accel 1\n');

            fprintf(fid,'recorder Node -file {%s} -timeSeries 1 -node 0 -dof 1 accel\n',filename_output_timeSeries);
            fprintf(fid,'recorder Node -file {%s} -time -nodeRange 1 %i -dof 1 disp \n',filename_output_def,obj.nStories);
            fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i force \n',filename_output_force,obj.nStories);
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
            fprintf(fid,'while { $currentTime < %g } {\n',tend);
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
            [status, result] = obj.runOpenSees(filename_input);
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
            temp = dlmread(filename_output_timeSeries);
            results.groundMotion = temp;
            temp = dlmread(filename_output_def);
            results.time = temp(:,1);
            results.totalDrift = temp(:,2:end);
            temp = dlmread(filename_output_force);
            results.storyShear = temp(:,2:2:end);

            % Computed Results
            storyDrift = results.totalDrift;
            storyDrift(:,2:end) = storyDrift(:,2:end)-storyDrift(:,1:(end-1));
            results.storyDrift = storyDrift;
            results.roofDrift = results.totalDrift(:,end);
            results.baseShear = results.storyShear(:,1);

            % Clean Folder
            if obj.deleteFilesAfterAnalysis
                delete(filename_input,filename_output_timeSeries,...
                    filename_output_def,filename_output_force);
            end
        end %function:responseHistory

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
            %% Story spring design

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
                spring(i).K0       = designStiffness(i);            % elastic stiffness
                spring(i).as       = springGivens.as;               % strain hardening ratio
                spring(i).Lambda_S = springGivens.Lambda_S;         % Cyclic deterioration parameter - strength
                spring(i).Lambda_K = springGivens.Lambda_K;         % Cyclic deterioration parameter - stiffness
                spring(i).c_S      = springGivens.c_S;              % rate of deterioration - strength
                spring(i).c_K      = springGivens.c_K;              % rate of deterioration - stiffness
                spring(i).Res      = springGivens.Res;              % residual strength ratio
                spring(i).D        = springGivens.D;                % rate of cyclic deterioration
                spring(i).nFactor  = springGivens.nFactor;          % elastic stiffness amplification factor

                spring(i).V_c = designStrength(i);                  % strength at capping
                spring(i).V_y = springGivens.C_yc*spring(i).V_c;    % effective yield strength

                spring(i).defl_y  = spring(i).V_y./spring(i).K0;                                                    % deflection at yield
                spring(i).defl_p  = (spring(i).V_c-spring(i).V_y)./(spring(i).as*spring(i).K0);                     % pre-capping deflection
                % spring(i).defl_pc = springGivens.C_pcp*spring(i).defl_p;                                            % post-capping deflection
                spring(i).defl_pc = spring(i).V_c/(springGivens.ad*spring(i).K0);
                % spring(i).defl_pc = spring(i).V_c*tand(springGivens.theta_pc);
                spring(i).defl_u  = springGivens.C_upc*(spring(i).defl_y + spring(i).defl_p + spring(i).defl_pc);   % ultimate deflection capacity

                spring(i).definition = bilinearMaterialDefinition(i,spring(i));
            end

        end %function:springDesign

        function animateResponseHistory(obj,results,dt)
            %% ANIMATERESPONSEHISTORY Animate a given response history
            %
            %   ANIMATERESPONSEHISTORY(obj,results) is the default usage.
            %
            %   ANIMATERESPONSEHISTORY(obj,results,dt) allows for overriding the
            %       timestep used. This option must be used to allow for
            %       pushover analyses, as those contain no time results and time
            %       is arbitrary with them anyways.
            %

            if nargin == 2
                if isfield(results,'time')
                    dt = max(diff(results.time));
                else
                    error('Undefined timestep.')
                end
            end

            cumHeights = cumsum(obj.storyHeight);

            figure

            grid on
            grid minor

            xMax = max(max(abs(results.totalDrift)));
            yMax = sum(obj.storyHeight) + obj.storyHeight(1);

            axis([-xMax xMax 0 yMax]);
            ax = gca;
            ax.YTick = cumHeights;

            ylabel(sprintf('Height (%s)',obj.units.length))
            xlabel(sprintf('Drift (%s)',obj.units.length))

            h = animatedline([0 results.totalDrift(1,:)],[0 cumHeights],'Marker','*');
            t = text(2/3*xMax,1/5*yMax,'Time: 0.0s');

            for i = 1:size(results.totalDrift,1)
                displayText = sprintf('Time: %4.1fs',i*dt);
                delete(t)
                clearpoints(h)
                addpoints(h,[0 results.totalDrift(i,:)],[0 cumHeights])
                t = text(ax,2/3*xMax,1/5*yMax,displayText);
                drawnow
            end

        end
    end %methods
end %classdef:mdofShearBuilding2d


function tf = isvectorsize(v,n)
    s = size(v);
    tf = isequal(s,[n 1]) || isequal(s,[1 n]);
end %function:isvectorsize

function writeFunction_updateRayleighDamping(fid)
    fprintf(fid,'proc updateRayleighDamping { modeA ratioA modeB ratioB } {\n');
    fprintf(fid,'    # ###################################################################\n');
    fprintf(fid,'    # updateRayleighDamping $modeA $ratioA $modeB $ratioB\n');
    fprintf(fid,'    # ###################################################################\n');
    fprintf(fid,'    # Runs an eigenvalue analysis and set proportional damping based on\n');
    fprintf(fid,'    # the current state of the structure\n');
    fprintf(fid,'    #\n');
    fprintf(fid,'    # Input Parameters:\n');
    fprintf(fid,'    # modeA, modeB - modes that will have prescribed damping ratios\n');
    fprintf(fid,'    # ratioA, ratioB - damping ratios prescribed at the specified modes\n');
    fprintf(fid,'\n');
    fprintf(fid,'    # Get natural frequencies at the desired modes\n');
    fprintf(fid,'    if { $modeA > $modeB } {\n');
    fprintf(fid,'        set maxMode $modeA\n');
    fprintf(fid,'    } else {\n');
    fprintf(fid,'        set maxMode $modeB\n');
    fprintf(fid,'    }\n');
    fprintf(fid,'\n');
    fprintf(fid,'    set eigs    [eigen -fullGenLapack $maxMode]\n');
    fprintf(fid,'    set freqA   [expr sqrt([lindex $eigs [expr $modeA-1]])]\n');
    fprintf(fid,'    set freqB   [expr sqrt([lindex $eigs [expr $modeB-1]])]\n');
    fprintf(fid,'\n');
    fprintf(fid,'    # Compute the damping factors\n');
    fprintf(fid,'    set tempVal [expr 2.0/($freqA*$freqA-$freqB*$freqB)]\n');
    fprintf(fid,'    set aM      [expr $tempVal*$freqA*$freqB*($ratioB*$freqA-$ratioA*$freqB)]\n');
    fprintf(fid,'    set aK      [expr $tempVal*($ratioA*$freqA-$ratioB*$freqB)]\n');
    fprintf(fid,'\n');
    fprintf(fid,'    # Set the damping\n');
    fprintf(fid,'    rayleigh $aM 0.0 0.0 $aK\n');
    fprintf(fid,'}\n');
end %function:writeFunction_updateRayleighDamping
