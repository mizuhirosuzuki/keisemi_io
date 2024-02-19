function Estimation_PS_bootstrap(
    FakeData,
    global_param
)

    EstimatedCCP1Mat, EstimatedCCP2Mat = estimateCCPMat(
        FakeData, global_param
    );

    EstimatedTransition = estimateTransition(FakeData);

    result = optimize(
        x -> calculatePSDObjective(
            x,
            EstimatedCCP1Mat, 
            EstimatedCCP2Mat, 
            EstimatedTransition, 
            global_param
            ),
        [0.3, 0.2, -0.27, 0.45, -2.1],
        Optim.Options(show_trace = false)
    );

    return [
        EstimatedCCP1Mat, 
        EstimatedCCP2Mat,
        EstimatedTransition,
        result.minimizer,
        result.minimum
    ]

end
