function Bootstrap_PS_forward(
    FakeData,
    global_param,
    EVrandom,
    UNIrandom,
    InitialState,
    NumSimMarkets,
    NumSimulations,
    NumSimPeriods
)

    EstimatedCCP1Mat, EstimatedCCP2Mat = estimateCCPMat(
        FakeData, global_param
    );

    EstimatedTransition = estimateTransition(FakeData);

    output = VSigmaGeneration(
        EstimatedCCP1Mat[:, 2],
        EstimatedCCP2Mat[:, 2],
        EstimatedTransition,
        EVrandom,
        UNIrandom,
        InitialState,
        NumSimMarkets,
        NumSimulations,
        NumSimPeriods,
        global_param
    );

    W1star = output[1];
    W2star = output[2];

    initial = [0.3, 0.2, -0.27, 0.45, -2.1];

    result = optimize(
        x -> Estimation_forward_PSD(
            x,
            W1star,
            W2star,
            EstimatedTransition,
            EstimatedCCP1Mat,
            EstimatedCCP2Mat,
            global_param
        ),
        initial,
        Optim.Options(show_trace = false)
    )

    return [
        result,
        EstimatedCCP1Mat,
        EstimatedCCP2Mat,
        EstimatedTransition
    ]
end
