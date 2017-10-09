function a = tern(cond, t, f)
%TERN  Ternary operator for variable assignment.
    if cond
        a = t;
    else
        a = f;
    end
end
