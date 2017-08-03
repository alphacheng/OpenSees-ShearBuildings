function tf = isFloatInt(x)
% ISFLOATINT Evaluate whether the input is an "integer", even if floating point
%
% Returns true if input is finite number with no trailing pieces:
% >> x = 3;
% >> class(x)
% ans =
%   'double'
% >> isFloatInt(x)
% ans =
%   logical
%    1
%
    tf = isnumeric(x) && isfinite(x) && (x == floor(x));
end
