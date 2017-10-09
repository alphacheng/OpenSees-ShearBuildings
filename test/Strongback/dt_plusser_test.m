dt = zeros(44,1);
spacing = zeros(44,1);
time = cell(44,1);
accel = cell(44,1);
for index = 1:44
    dt(index) = ground_motions(index).dt;
    if dt(index) < 0.01
        spacing(index) = floor(0.01/dt(index));
        dt(index) = 0.01;
    else
        spacing(index) = 1;
    end
    time{index} = ground_motions(index).time(1:spacing(index):end);
    accel{index} = ground_motions(index).normalized_acceleration(1:spacing(index):end);
end
