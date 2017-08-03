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

assert(isnumeric(F) && isvector(F) && (length(F) == obj.nStories),...
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
filenames.input              = obj.scratchFile('mdofShearTrussBuilding_input.tcl');
filenames.output_timeSeries  = obj.scratchFile('mdofShearTrussBuilding_timeSeries.out');
filenames.output_def_x       = obj.scratchFile('mdofShearTrussBuilding_disp_x.out');
filenames.output_def_y       = obj.scratchFile('mdofShearTrussBuilding_disp_y.out');
filenames.output_force_story = obj.scratchFile('mdofShearTrussBuilding_force_s.out');
filenames.output_force_truss = obj.scratchFile('mdofShearTrussBuilding_force_t.out');

% Create .tcl file
fid = fopen(filenames.input,'w');

obj.constructBuilding(fid)
if obj.includeExplicitPDelta
    obj.applyGravityLoads(fid)
end

fprintf(fid,'timeSeries Linear 1\n');
fprintf(fid,'pattern Plain 1 1 {\n');
for i = 1:obj.nStories
    if F(i) ~= 0
        fprintf(fid,'    load %i %g 0\n',i*10,F(i));
    end
end
fprintf(fid,'} \n');

fprintf(fid,'recorder Node -file {%s} -time -timeSeries 1 -node 0 -dof 1 accel\n',filenames.output_timeSeries);
fprintf(fid,'recorder Node -file {%s} -nodeRange 10 %i -dof 1 disp \n',filenames.output_def_x,obj.nStories*10);
fprintf(fid,'recorder Node -file {%s} -nodeRange 10 %i -dof 2 disp \n',filenames.output_def_y,obj.nStories*10);
fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i force \n',filenames.output_force_story,obj.nStories);
fprintf(fid,'recorder Element -file {%s} -eleRange 10 %i force \n',filenames.output_force_truss,obj.nStories*10);
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
fprintf(fid,'integrator DisplacementControl %i 1 %g\n',controlStory1(obj),obj.optionsPushover.stepSize);
fprintf(fid,'analysis Static \n');

switch lower(type)
    case 'targetdrift'
        fprintf(fid,'algorithm Newton\n');
        fprintf(fid,'set ok [analyze %i]\n',ceil(targetDrift/obj.optionsPushover.stepSize));
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
        fprintf(fid,'    set currentLoad [lindex [eleResponse 1 force] 2]\n');
        % fprintf(fid,'    set currentLoad [getTime]\n');
        fprintf(fid,'    if { $currentLoad > $maxLoad } {\n');
        fprintf(fid,'        set maxLoad $currentLoad\n');
        fprintf(fid,'    }\n');
        fprintf(fid,'    if { [nodeDisp %i 1] > %g } {\n',controlStory1(obj),obj.optionsPushover.maxDrift);
        fprintf(fid,'        exit 3\n');
        fprintf(fid,'    }\n');
        fprintf(fid,'}\n');
    otherwise
        error('Unknown analysis type: %s',type);
end

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
    case 3
        results.exitStatus = 'Peak Drift Reached';
    otherwise
        fprintf('%s\n',result);
        error('Analysis Failed in Unknown Manner (exit code: %i)',status);
end

% Read Results
temp = dlmread(filenames.output_timeSeries);
time = temp(:,1);
results.displacement_x = dlmread(filenames.output_def_x);
results.displacement_y = dlmread(filenames.output_def_y);

temp = dlmread(filenames.output_force_story);
results.storyShear = -temp(:,1:4:end);

temp = dlmread(filenames.output_force_truss);
results.trussForce_x = temp(:,1:4:end);
results.trussForce_y = temp(:,2:4:end);

% Computed Results
storyDrift = results.displacement_x;
storyDrift(:,2:end) = storyDrift(:,2:end)-storyDrift(:,1:(end-1));
results.storyDrift = storyDrift;
results.roofDrift = results.displacement_x(:,end);
results.appliedStoryForce = time*F;
results.baseShear = sum(results.appliedStoryForce,2);

% Clean Folder
if obj.deleteFilesAfterAnalysis
    fields = fieldnames(filenames);
    for i = 1:length(fields)
        delete(filenames.(fields{i}));
    end
end

function s = controlStory1(obj)
    if isequal(obj.optionsPushover.controlStory,'roof')
        s = obj.nStories*10;
    else
        s = obj.optionsPushover.controlStory*10;
    end
end

end %function:pushover
