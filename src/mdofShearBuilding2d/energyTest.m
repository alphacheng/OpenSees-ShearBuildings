function E_EQ = energyTest(obj,results)
% ENERGYTEST

time = results.time;

M = obj.storyMass(:);

relativeVeloc = results.totalVeloc;

u_dot = relativeVeloc;
u_ddot_eq = results.groundMotion;

expression = u_dot*M.*u_ddot_eq;

E_EQ = -cumtrapz(time,expression);

end
