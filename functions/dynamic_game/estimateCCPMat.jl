function estimateCCPMat(data, global_param)

    EstimatedCCP1 = zeros(8, 1);
    EstimatedCCP2 = zeros(8, 1);

    for s in 1:8
        SubData = data[data[:, 3] .== s, :];
        EstimatedCCP1[s] = mean(SubData[:, 7] .== 0);
        EstimatedCCP2[s] = mean(SubData[:, 8] .== 0);
    end

    EstimatedCCP1Mat = (
        hcat(
            1 .- EstimatedCCP1,
            EstimatedCCP1,
            1 .- EstimatedCCP1,
        ) .* 
        global_param.possibleActionArray[:, :, 1]
        );
    EstimatedCCP2Mat = (
        hcat(
            1 .- EstimatedCCP2,
            EstimatedCCP2,
            1 .- EstimatedCCP2,
        ) .* 
        global_param.possibleActionArray[:, :, 2]
        );

    return [EstimatedCCP1Mat, EstimatedCCP2Mat]
end
