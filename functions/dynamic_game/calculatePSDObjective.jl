function calculatePSDObjective(
    param,
    EstimatedCCP1Mat,
    EstimatedCCP2Mat,
    TransitionProb,
    global_param
)

    parameterMat = hcat(
        [param[1]; param[3:4]; [0]; param[5]],
        [param[2]; param[3:4]; [0]; param[5]], 
    );

    output = updateCCP(
        parameterMat, 
        EstimatedCCP1Mat, 
        EstimatedCCP2Mat, 
        TransitionProb, 
        global_param
    );

    return (
        sum((output[1][:, 2] - EstimatedCCP1Mat[:, 2]).^2) + 
        sum((output[2][:, 2] - EstimatedCCP2Mat[:, 2]).^2)
    )

end;
