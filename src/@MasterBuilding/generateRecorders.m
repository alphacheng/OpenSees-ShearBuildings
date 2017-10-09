function recorders = generateRecorders(obj, analysis, filenames)
% GENERATERECORDERS

recorders = struct;
roof = obj.nStories;

recorders.timeSeries = sprintf('Node -file {%s} -time -timeSeries 1 -node %i -dof 1 accel', filenames.output_timeSeries, obj.storyNodes(0));
recorders.disp_x     = sprintf('Node -file {%s} -nodeRange %i %i -dof 1 disp', filenames.output_disp_x, obj.storyNodes(1), obj.storyNodes(roof));
recorders.disp_y     = sprintf('Node -file {%s} -nodeRange %i %i -dof 2 disp', filenames.output_disp_y, obj.storyNodes(1), obj.storyNodes(roof));
recorders.vel_x      = sprintf('Node -file {%s} -nodeRange %i %i -dof 1 vel',  filenames.output_vel_x,  obj.storyNodes(1), obj.storyNodes(roof));
recorders.vel_y      = sprintf('Node -file {%s} -nodeRange %i %i -dof 2 vel',  filenames.output_vel_y,  obj.storyNodes(1), obj.storyNodes(roof));
recorders.force_s    = sprintf('Element -file {%s} -eleRange %i %i -dof 1 force', filenames.output_force_story, obj.springTags(1), obj.springTags(roof));
recorders.force_t    = sprintf('Element -file {%s} -eleRange %i %i -dof 1 2 force', filenames.output_force_truss, obj.trussTags(1), obj.trussTags(roof));
recorders.force_b    = sprintf('Element -file {%s} -eleRange %i %i -dof 1 2 force', filenames.output_force_sback, obj.sbackTags(1), obj.)



fprintf(fid,'recorder Node -file {%s} -time -timeSeries 1 -node %i -dof 1 accel\n', filenames.output_timeSeries, obj.storyNodes(0));
fprintf(fid,'recorder Node -file {%s} -nodeRange %i %i -dof 1 disp \n', filenames.output_def_x, obj.storyNodes(1), obj.storyNodes(obj.nStories));
fprintf(fid,'recorder Node -file {%s} -nodeRange %i %i -dof 2 disp \n', filenames.output_def_y, obj.storyNodes(1), obj.storyNodes(obj.nStories));
fprintf(fid,'recorder Element -file {%s} -eleRange 1 %i -dof 1 force \n', filenames.output_force_story, obj.nStories);
fprintf(fid,'recorder Element -file {%s} -eleRange %i %i -dof 1 2 force \n',filenames.output_force_truss, obj.trussTags(1), obj.trussTags(obj.nStories));
fprintf(fid,'recorder Element -file {%s} -eleRange 21 %i -dof 1 2 force \n',filenames.output_force_sback,obj.nStories+20);


recorders = obj.generateRecorders('pushover', filenames);
switch class(obj)
case 'mdofShearBuilding2d_new'
    fprintf(fid,'recorder %s\n', recorders.timeSeries);
    fprintf(fid,'recorder %s\n', recorders.disp_x);
    fprintf(fid,'recorder %s\n', recorders.vel_x);
    fprintf(fid,'recorder %s\n', recorders.force_s);
case 'mdofShearTrussBuilding'
    fprintf(fid,'recorder %s\n', recorders.timeSeries);
    fprintf(fid,'recorder %s\n', recorders.disp_x);
    fprintf(fid,'recorder %s\n', recorders.disp_y);
    fprintf(fid,'recorder %s\n', recorders.vel_x);
    fprintf(fid,'recorder %s\n', recorders.vel_y);
    fprintf(fid,'recorder %s\n', recorders.force_s);
    fprintf(fid,'recorder %s\n', recorders.force_t);
case 'Strongback'
    fprintf(fid,'recorder %s\n', recorders.timeSeries);
    fprintf(fid,'recorder %s\n', recorders.disp_x);
    fprintf(fid,'recorder %s\n', recorders.disp_y);
    fprintf(fid,'recorder %s\n', recorders.vel_x);
    fprintf(fid,'recorder %s\n', recorders.vel_y);
    fprintf(fid,'recorder %s\n', recorders.force_s);
    fprintf(fid,'recorder %s\n', recorders.force_t);
    fprintf(fid,'recorder %s\n', recorders.force_b);
end
