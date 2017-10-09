classdef MasterBuilding < OpenSeesAnalysis

properties

%------------------------------ General settings ------------------------------%

g                               % Acceleration due to gravity
units                           % Units used for labels
verbose = true;                 % Print updates on progress to MATLAB console
pathOf                          % Paths to OpenSees, external functions, etc.


%----------------------------- Building definition ----------------------------%

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


%------------------------------ Analysis options ------------------------------%

% optionsGravityLoads - class containing settings for applying gravity loads
optionsGravityLoads = OpenSees.GravityLoadOptions();

% optionsPushover - class containing settings for pushover analysis
optionsPushover = OpenSees.PushoverOptions();

% optionsResponseHistory - class containing settings for response history analysis
optionsResponseHistory = OpenSees.ResponseHistoryOptions();

% optionsIDA - class containing settings for incremental dynamic analysis
optionsIDA = OpenSees.IncrementalDynamicAnalysisOptions();


%---------------------- Equivalent Lateral Force options ----------------------%

seismicDesignCategory       % Seismic design category (ASCE 7-10 Section 11.6)
impFactor                   % Seismic importance factor (ASCE 7-10 Section 1.5)
respModCoeff                % Response modification coefficient (ASCE 7-10 Section 12.2)
deflAmplFact                % Deflection amplification factor (ASCE 7-10 Section 12.2)
overstrengthFactor          % Overstrength factor (ASCE 7-10 Section 12.2)

end

methods
function obj = MasterBuilding(nStories)
    obj.nStories = nStories;
end

%############################# Set et get methods #############################%
%#ok<*MCSUP> -- Code analyzer complains about accessing other properties in
% these set methods; this isn't a problem because the property being accessed
% (nStories) is set by the constructor.

function set.storyMass(obj, storyMass)
    assert(isnumeric(storyMass) && isvector(storyMass) && (length(storyMass) == obj.nStories),...
            'storyMass must be a numeric vector of length %i', obj.nStories)
    obj.storyMass = storyMass;
end
function set.storyHeight(obj, storyHeight)
    assert(isnumeric(storyHeight) && isvector(storyHeight) && (length(storyHeight) == obj.nStories),...
            'storyHeight must be a numeric vector of length %i', obj.nStories)
    obj.storyHeight = storyHeight;
end
function set.storySpringDefinition(obj, ssdef)
    assert(isvector(ssdef) && (length(ssdef) == obj.nStories) && all(cellfun(@ischar, ssdef)),...
            'storySpringDefinition must be a cell vector, length %i, of char vectors', obj.nStories)
    obj.storySpringDefinition = ssdef;
end
function set.storyTrussDefinition(obj, stdef)
    assert(isvector(stdef) && (length(stdef) == obj.nStories) && all(cellfun(@ischar, stdef)),...
            'storyTrussDefinition must be a cell vector, length %i, of char vectors', obj.nStories)
    obj.storyTrussDefinition = stdef;
end

%################################## Analyses ##################################%

function [vals, vec] = eigenvalues(obj, n, a)
%% EIGENVALUES Eigenvalue analysis of system.
%
%   vals = EIGENVALUES(obj, n) returns the first n eigenvalues of obj in the
%       vector vals. Note that the eigenvalues are equal to the square of the
%       circular natural frequencies, not the frequencies themselves. If
%       unspecified, n defaults to the number of stories in obj.
%
%   [vals, vec] = EIGENVALUES(obj, n, a) returns the first n eigenvalues of obj
%       in the vector vals and the eigenvector of the a-th mode of all nodes of
%       obj in the vector vec. If unspecified, a defaults to 1.
%
    if nargin < 3
        a = 1;
    end
    if nargin < 2
        n = obj.nStories;
    end

    if n == obj.nStories
        solver = '-fullGenLapack ';
    elseif n < obj.nStories
        solver = '';
    else
        error('n must be less than or equal to the number of stories')
    end

    if a > n
        error('Insufficient eigenvalues (%i) to solve for requested eigenvector (%i)', n, a)
    end

    filenames = obj.generateFilenames('eigen');

    fid = fopen(filenames.input,'w');

    obj.constructBuilding(fid)
    if obj.includeExplicitPDelta
        obj.applyGravityLoads(fid)
    end

    fprintf(fid,'############################## Eigenvalue Analysis #############################\n');
    fprintf(fid,'set eigs [eigen %s%i]\n', solver, n);
    fprintf(fid,'set eigenvalues $eigs\n');
    fprintf(fid,'set vecs {}\n');
    fprintf(fid,'set vecfid [open %s w+]\n',obj.path_for_tcl(filenames.vecs));
    for i = 1:obj.nStories
        fprintf(fid,'lappend vecs [nodeEigenvector %i %i]', obj.storyNodes(i), a);
        fprintf(fid,'puts $vecfid [lindex $vecs %i 0]', i-1);
    end
    % fprintf(fid,'for {set i 1} {$i <= %i} {incr i} {\n',obj.nStories);
    % fprintf(fid,'  lappend vecs [nodeEigenvector [expr $i+10] %i]\n', a);
    % fprintf(fid,'  puts $vecfid [lindex $vecs [expr $i - 1] 0]\n');
    % fprintf(fid,'}\n');
    fprintf(fid,'set eigfid [open %s w+]\n',obj.path_for_tcl(filenames.vals));
    fprintf(fid,'puts $eigfid $eigs\n');
    fprintf(fid,'close $eigfid\n');
    fprintf(fid,'close $vecfid\n');
    fclose(fid);

    [~,~] = obj.runOpenSees(filenames.input);

    vals = dlmread(filenames.vals);
    if nargout == 2
        vec = dlmread(filenames.vecs);
    end

    if obj.deleteFilesAfterAnalysis
        delete(filenames.input,filenames.vals,filenames.vecs);
    end

end %function:eigenvalues

end

end
