classdef Strongback < OpenSeesAnalysis

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
    strongbackDefinition            % Struct containing `Area`, `Modulus`, and `Inertia` for the strongback

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
    optionsGravityLoads = struct('constraints' , OpenSees.ConstraintsOptions(), ...
                                 'test'        , OpenSees.TestOptions('tolerance',1e-6), ...
                                 'algorithm'   , OpenSees.AlgorithmOptions('KrylovNewton') ...
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
                             'constraints' , OpenSees.ConstraintsOptions('type','Plain'), ...
                             'test'        , OpenSees.TestOptions(), ...
                             'algorithm'   , OpenSees.AlgorithmOptions() ...
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
    optionsResponseHistory = struct('constraints' , OpenSees.ConstraintsOptions(), ...
                                    'test'        , OpenSees.TestOptions('tolerance',[1e-5,5e-5,1e-4]), ...
                                    'algorithm'   , OpenSees.AlgorithmOptions() ...
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

end %properties

methods

%###############################################################################
%% Constructor and set methods #################################################
%###############################################################################

function obj = Strongback(nStories)
    obj.nStories = nStories;
end

%###############################################################################
%% Shared model functions ######################################################
%###############################################################################

% fprintf(fid,'source [file join {%s} {updateRayleighDamping.tcl}]\n',obj.pathOf.tclfunctions);

function constructBuilding(obj,fid)
%% CONSTRUCTBUILDING Create the OpenSees model based on current properties
%
%    CONSTRUCTBUILDING(obj,fid) writes the OpenSees code that represents
%       the model to the file specified by fid.
%
    cumHeights = cumsum(obj.storyHeight);

    fprintf(fid,'# Units: %s, %s, %s\n\n',obj.units.force,obj.units.length,obj.units.time);
    fprintf(fid,'#################################### Model #####################################\n');
    fprintf(fid,'model BasicBuilder -ndm 2 -ndf 3 \n\n');
    fprintf(fid,'#----------------------------------- Nodes ------------------------------------#\n');
    fprintf(fid,'# Ground\n');
    fprintf(fid,'node 10 0 0\n');
    fprintf(fid,'node 20 0 0\n');
    fprintf(fid,'# Springs\n');
    for i = 1:obj.nStories
        fprintf(fid,'node %i 0 %g\n',i,cumHeights(i));
    end
    fprintf(fid,'# Trusses\n');
    for i = 1:obj.nStories
        fprintf(fid,'node %i 0 %g -mass %g %g %g\n',10+i,cumHeights(i),obj.storyMass(i),obj.storyMass(i),obj.storyMass(i));
    end
    fprintf(fid,'# Strongback\n');
    for i = 1:obj.nStories
        fprintf(fid,'node %i 0 %g\n',20+i,cumHeights(i));
    end
    fprintf(fid,'\n');
    fprintf(fid,'#-------------------------------- Constraints ---------------------------------#\n');
    fprintf(fid,'fix 10 1 1 1\n');
    fprintf(fid,'fix 20 1 1 0\n\n');
    fprintf(fid,'# Truss nodes need to be fixed in rotation\n');
    for i = 1:obj.nStories
        fprintf(fid,'fix %i 0 0 1\n',10+i);
    end
    fprintf(fid,'\n');
    fprintf(fid,'# Attach springs to trusses\n');
    for i = 1:obj.nStories
        fprintf(fid,'equalDOF %i %i 1 2 3\n',10+(i-1),i);
    end
    fprintf(fid,'\n');
    fprintf(fid,'# Attach trusses to strongback\n');
    for i = 1:obj.nStories
        fprintf(fid,'equalDOF %i %i 1\n',10+i,20+i);
    end
    fprintf(fid,'\n');
    fprintf(fid,'#--------------------------------- Materials ----------------------------------#\n');
    fprintf(fid,'# Springs\n');
    for i = 1:length(obj.storySpringDefinition)
        fprintf(fid,'%s\n',obj.storySpringDefinition{i});
    end
    fprintf(fid,'\n');
    fprintf(fid,'# Trusses\n');
    for i = 1:length(obj.storyTrussDefinition)
        fprintf(fid,'%s\n',obj.storyTrussDefinition{i});
    end
    fprintf(fid,'\n');
    fprintf(fid,'#---------------------------------- Elements ----------------------------------#\n');
    fprintf(fid,'# Springs\n');
    for i = 1:obj.nStories
        fprintf(fid,'element zeroLength %i %i %i -mat %i -dir 1 -doRayleigh 1\n',i,i,i+10,i);
    end
    fprintf(fid,'\n');
    fprintf(fid,'# Trusses\n');
    for i = 1:obj.nStories
        fprintf(fid,'element corotTruss %i %i %i 1 %i\n',i+10,(i-1)+10,i+10,i+10);
    end
    fprintf(fid,'\n');
    fprintf(fid,'# Strongback\n');
    fprintf(fid,'geomTransf Corotational 1\n');
    sbArea    = obj.strongbackDefinition.Area;
    sbModulus = obj.strongbackDefinition.Modulus;
    sbInertia = obj.strongbackDefinition.Inertia;
    for i = 1:obj.nStories
        fprintf(fid,'element elasticBeamColumn %i %i %i %g %g %g 1\n',20+i,20+(i-1),20+i,sbArea,sbModulus,sbInertia);
    end
    fprintf(fid,'\n');

end %function:constructBuilding

function applyGravityLoads(obj,fid)
%% APPLYGRAVITYLOADS Write gravity load commands to file

    fprintf(fid,'################################# Gravity Loads ################################\n');
    fprintf(fid,'pattern Plain 0 Linear {\n');
    for i = 1:obj.nStories
        fprintf(fid,'    load %i 0 -%g 0\n',i+10,obj.storyMass(i)*obj.g);
    end
    fprintf(fid,'}\n\n');
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
    fprintf(fid,'algorithm %s\n',obj.optionsGravityLoads.algorithm.type{1});
    fprintf(fid,'integrator LoadControl 0.1\n');
    fprintf(fid,'analysis Static\n\n');
    fprintf(fid,'analyze 10\n');
    fprintf(fid,'loadConst -time 0.0\n');
    fprintf(fid,'wipeAnalysis\n\n');
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

    filename_input = obj.scratchFile('Strongback_eigen_input.tcl');
    filename_vals  = obj.scratchFile('Strongback_eigen_vals.out');
    filename_vecs  = obj.scratchFile('Strongback_eigen_vecs.out');

    fid = fopen(filename_input,'w');

    obj.constructBuilding(fid)
    if obj.includeExplicitPDelta
        obj.applyGravityLoads(fid)
    end

    fprintf(fid,'############################## Eigenvalue Analysis #############################\n');
    fprintf(fid,'set eigs [eigen -fullGenLapack %i]\n',obj.nStories);
    fprintf(fid,'set eigenvalues $eigs\n');
    fprintf(fid,'set vecs {}\n');
    fprintf(fid,'set vecfid [open %s w+]\n',obj.path_for_tcl(filename_vecs));
    fprintf(fid,'for {set i 1} {$i <= %i} {incr i} {\n',obj.nStories);
    fprintf(fid,'  lappend vecs [nodeEigenvector [expr $i+10] 1]\n');
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

end %methods

end %classdef
