load '../ground_motions.mat'
dt_min = 0.01;

maxPosDiff = zeros(44,1);
maxNegDiff = zeros(44,1);
maxDiff = zeros(44,1);

for i = 1:44
    dt = ground_motions(i).dt;
    if dt < dt_min
        spacing = floor(dt_min/dt);
        dt = dt_min;
    else
        spacing = 1;
    end
    accel = ground_motions(i).normalized_acceleration(1:spacing:end);
    maxPosDiff(i) = max(ground_motions(i).normalized_acceleration) - max(accel);
    maxNegDiff(i) = min(ground_motions(i).normalized_acceleration) - min(accel);
    maxDiff(i) = max(abs(maxPosDiff(i)), abs(maxNegDiff(i)));
end

maximum_accel_diff = max(maxDiff);
