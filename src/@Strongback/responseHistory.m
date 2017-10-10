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
filenames.input              = obj.scratchFile(sprintf('%s_responseHistory_input_%s_%i.tcl'      ,class(obj),gmID,indexNum));
filenames.output_timeSeries  = obj.scratchFile(sprintf('%s_responseHistory_timeSeries_%s_%i.out' ,class(obj),gmID,indexNum));
filenames.output_disp_truss  = obj.scratchFile(sprintf('%s_responseHistory_disp_truss_%s_%i.out' ,class(obj),gmID,indexNum));
filenames.output_disp_sback  = obj.scratchFile(sprintf('%s_responseHistory_disp_sback_%s_%i.out' ,class(obj),gmID,indexNum));
filenames.output_vel_x       = obj.scratchFile(sprintf('%s_responseHistory_vel_x_%s_%i.out'      ,class(obj),gmID,indexNum));
filenames.output_vel_y       = obj.scratchFile(sprintf('%s_responseHistory_vel_y_%s_%i.out'      ,class(obj),gmID,indexNum));
filenames.output_force_story = obj.scratchFile(sprintf('%s_responseHistory_force_s_%s_%i.out'    ,class(obj),gmID,indexNum));
filenames.output_force_truss = obj.scratchFile(sprintf('%s_responseHistory_force_t_%s_%i.out'    ,class(obj),gmID,indexNum));
filenames.output_force_sback = obj.scratchFile(sprintf('%s_responseHistory_force_b_%s_%i.out'    ,class(obj),gmID,indexNum));
filenames.output_reaction    = obj.scratchFile(sprintf('%s_responseHistory_reaction_%s_%i.out'   ,class(obj),gmID,indexNum));

% Create .tcl file
fid = fopen(filenames.input,'w');

obj.constructBuilding(fid)
if obj.includeExplicitPDelta
    obj.applyGravityLoads(fid)
end

fprintf(fid,'############################### Response history ###############################\n');
fprintf(fid,'timeSeries Path 1 -dt %g -filePath {%s} -factor %g\n',dt,gmFile,SF);
fprintf(fid,'pattern UniformExcitation 1 1 -accel 1\n\n');

fprintf(fid,'#---------------------------------- Recorders ---------------------------------#\n');
fprintf(fid,'recorder Node -file {%s} -time -timeSeries 1 -node 10 -dof 1 accel\n',filenames.output_timeSeries);
fprintf(fid,'recorder Node -file {%s} -nodeRange 11 %i -dof 1 2 disp \n',filenames.output_disp_truss,obj.nStories+10);
fprintf(fid,'recorder Node -file {%s} -nodeRange 21 %i -dof 2 disp \n',filenames.output_disp_sback,obj.nStories+20);
fprintf(fid,'recorder Node -file {%s} -nodeRange 11 %i -dof 1 vel \n',filenames.output_vel_x,obj.nStories+10);
fprintf(fid,'recorder Node -file {%s} -nodeRange 11 %i -dof 2 vel \n',filenames.output_vel_y,obj.nStories+10);
fprintf(fid,'recorder Node -file {%s} -node 10 20 -dof 1 2 3 reaction \n',filenames.output_reaction);
fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i -dof 1 force \n',filenames.output_force_story,obj.nStories);
fprintf(fid,'recorder Element -file {%s} -eleRange 11 %i -dof 1 2 force \n',filenames.output_force_truss,obj.nStories+10);
fprintf(fid,'recorder Element -file {%s} -eleRange 21 %i localForce \n',filenames.output_force_sback,obj.nStories+20);
fprintf(fid,'record \n\n');

fprintf(fid,'#---------------------------------- Analysis ----------------------------------#\n');
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
fprintf(fid,'    algorithm %s\n',obj.optionsResponseHistory.algorithm.type{1});
fprintf(fid,'    test %s\n',testArgs{1});
fprintf(fid,'    set ok [analyze 1 %g]\n',dt);
for i = 1:length(testArgs)
    if i == 1, k = 2; else, k = 1; end
    for j = k:length(obj.optionsResponseHistory.algorithm.type)
        fprintf(fid,'    if { $ok != 0 } {\n');
        fprintf(fid,'        algorithm %s\n',obj.optionsResponseHistory.algorithm.type{j});
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

%-------------------------------- Read Results --------------------------------%
temp = dlmread(filenames.output_timeSeries);
results.time = temp(:,1);
results.groundMotion = temp(:,2);

temp = dlmread(filenames.output_disp_truss);
results.displacement_x = temp(:,1:2:end);
results.displacement_y = temp(:,2:2:end);

results.strongbackDispY = dlmread(filenames.output_disp_sback);

results.velocity_x     = dlmread(filenames.output_vel_x);
results.velocity_y     = dlmread(filenames.output_vel_y);

temp = dlmread(filenames.output_force_story);
results.storyShear = -temp;

temp = dlmread(filenames.output_force_truss);
results.trussForce_x = temp(:,1:2:end);
results.trussForce_y = temp(:,2:2:end);

temp = dlmread(filenames.output_force_sback);
results.strongbackAxial           = temp(:,1:6:end);
results.strongbackAxial(:,end+1)  = temp(:,end-2);
results.strongbackShear           = temp(:,2:6:end);
results.strongbackShear(:,end+1)  = temp(:,end-1);
results.strongbackMoment          = temp(:,3:6:end);
results.strongbackMoment(:,end+1) = temp(:,end);

reaction = dlmread(filenames.output_reaction);

%------------------------------ Computed Results ------------------------------%
storyDrift = results.displacement_x;
storyDrift(:,2:end) = storyDrift(:,2:end)-storyDrift(:,1:(end-1));
results.storyDrift = storyDrift;
results.roofDrift = results.displacement_x(:,end);
results.baseShear = results.storyShear(:,1) + reaction(:,1) + reaction(:,4);

% Clean Folder
if obj.deleteFilesAfterAnalysis
    fields = fieldnames(filenames);
    for i = 1:length(fields)
        delete(filenames.(fields{i}));
    end
end
end %function:responseHistory
