function estimateTransition(data)
    EstimatedTransition = zeros(2, 2);

    SubData = hcat(data, [[0]; data[1:(end - 1), 4]]);
    SubData = SubData[SubData[:, 2] .!= 1, :];

    for z in 1:2
        SubDataZ = SubData[SubData[:, 4] .== z, :];
        EstimatedTransition[z, z] = mean(SubDataZ[:, 9] .== z);
        EstimatedTransition[z, 3 - z] = 1 - EstimatedTransition[z, z];
    end

    SubData = hcat(data, [[0]; data[1:(end - 1), 4]]);
    SubData = SubData[SubData[:, 2] .!= 1, :];

    for z in 1:2
        SubDataZ = SubData[SubData[:, 4] .== z, :];
        EstimatedTransition[z, z] = mean(SubDataZ[:, 9] .== z);
        EstimatedTransition[z, 3 - z] = 1 - EstimatedTransition[z, z];
    end

    return EstimatedTransition
end
