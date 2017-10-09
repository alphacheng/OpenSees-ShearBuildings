classdef GravityLoadOptions < handle

properties
    constraints = OpenSees.ConstraintsOptions();
    test = OpenSees.TestOptions('tolerance',1e-6);
    algorithm = OpenSees.AlgorithmOptions('KrylovNewton');
end

methods

function set.constraints(obj,constraints)
    assert(isa(constraints,'ConstraintsOptions'), 'constraints must be a ConstraintsOptions object')
    obj.constraints = constraints;
end
function set.test(obj,test)
    assert(isa(test,'TestOptions'), 'test must be a TestOptions object')
    obj.test = test;
end
function set.algorithm(obj,algorithm)
    assert(isa(algorithm,'AlgorithmOptions'), 'algorithm must be a AlgorithmOptions object')
    obj.algorithm = algorithm;
end

end

end
