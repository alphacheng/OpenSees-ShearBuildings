function filenames = generateFilenames(obj, analysis, gmID, index)
% GENERATEFILENAMES  Generate a struct of consistently formatted file names

filenames = struct;
classname = class(obj);
prefix = sprintf('%s_%s_',classname,analysis);
if nargin == 4
    suffix = sprintf('_%s_%i.out', gmID, index);
else
    suffix = '.out';
end

filenames.input = 'input';
switch analysis
case 'eigen'
    filenames.vals = 'vals';
    filenames.vecs = 'vecs';
otherwise
    switch analysis
    case 'responseHistory'
        filenames.output_vel_x = 'vel_x';
        filenames.output_vel_y = 'vel_y';
    end
    filenames.output_timeSeries  = 'timeSeries';
    filenames.output_disp_x      = 'disp_x';
    filenames.output_disp_y      = 'disp_y';
    filenames.output_force_story = 'force_s';
    filenames.output_force_truss = 'force_t';
    filenames.output_force_sback = 'force_b';
end

filenames = structfun(@(x)sprintf('%s%s%s', prefix, x, suffix), filenames, 'UniformOutput', false);
filenames.input = strrep(filenames.input, '.out', '.tcl');
filenames = structfun(@(x)scratchFile(obj,x), filenames, 'UniformOutput', false);
filenames = orderfields(filenames);

end
